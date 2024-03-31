import gradio as gr
import os
from typing import Tuple, Dict, List
import pandas as pd
import geopandas as gpd
import json
import platform
import asyncio
import re
import glob

from osgeo import gdal, ogr

gdal.UseExceptions()

SYSTEM = platform.system()

with open('defaults_and_docs.json', 'r', encoding='utf-8') as f:
    d = json.loads(f.read())[0]
    docs: Dict[str:str] = d['docs']

def _format_files(file_path: str) -> str:
    """
    Function added so that windows paths that pad with quotation marks can be valid
    """
    if any(file_path.endswith(s) and file_path.startswith(s) for s in ['"',"'"]):
        return file_path[1:-1]
    return file_path

def _write(f, Card, Argument = '') -> None:
    if Argument:
        f.write(f"{Card}\t{Argument}\n")
    else:
        f.write(f"{Card}\n")

def _warn_DNE( card: str, value: str) -> None:
    if value:
        if not os.path.exists(value):
            gr.Warning(f"The file {value} for {card} does not exist: ")

def _check_type( card: str, value: str, types: List[str]) -> None:
    if not any(value.endswith(t) for t in types):
        gr.Warning(f"{card} is {os.path.basename(value)}, which does not have a valid file type ({','.join(types)})")

def omit_outliers_change(value: str) -> Tuple[gr.Radio, gr.Column, gr.Number]:
    if value == 'None':
        return gr.Radio(info='None: No outliers will be removed'), gr.Column(visible=False), gr.Number(visible=False)
    elif value == 'Flood Bad Cells':
        return gr.Radio(info=docs['Flood_BadCells']), gr.Column(visible=False), gr.Number(visible=False)
    elif value == 'Use AutoRoute Depths':
        return gr.Radio(info=docs['FloodSpreader_Use_AR_Depths']), gr.Column(visible=False), gr.Number(visible=False)
    elif value == 'Smooth Water Surface Elevation':
        return  gr.Radio(info=docs['smooth_wse']), gr.Column(visible=True), gr.Number(visible=False)
    elif value == 'Use AutoRoute Depths (StDev)':
        return gr.Radio(info=docs['FloodSpreader_Use_AR_Depths_StDev']), gr.Column(visible=False), gr.Number(visible=False)
    else:
        return gr.Radio(info=docs['FloodSpreader_SpecifyDepth']), gr.Column(visible=False), gr.Number(visible=True)

def bathy_changes(value: str) -> Tuple[gr.Slider, gr.Slider]:
    if value == 'Trapezoidal':
        return gr.Slider(visible=True), gr.Slider(visible=False)
    if value in ['Left Bank Quadratic', 'Right Bank Quadratic', 'Double Quadratic']:
        return gr.Slider(visible=True), gr.Slider(visible=True)
    return gr.Slider(visible=False), gr.Slider(visible=False)

def show_mans_n(lu, mannings) -> gr.Number:
    """
    Let mannings n input box be interactive if these two are both not specified
    """
    if lu and mannings:
        return gr.Number(interactive=False)
    return gr.Number(interactive=True)

def update_flow_params(flow_file) -> Tuple[gr.Dropdown, gr.Dropdown, gr.Dropdown, gr.Dropdown]:
    flow_file = _format_files(flow_file)
    if not os.path.exists(flow_file):
        return None, None, None, None
    try:
        cols = sorted(list(pd.read_csv(flow_file, delimiter='\t').columns))
        if len(cols) <= 1:
            try:
                cols = sorted(list(pd.read_csv(flow_file, delimiter=' ').columns))
            except:
                return None, None, None, None
    except:
        try:
            cols = sorted(list(pd.read_csv(flow_file, delimiter=' ').columns))
        except:
            return None, None, None, None

    return gr.Dropdown(choices=cols), gr.Dropdown(choices=cols), gr.Dropdown(choices=cols), gr.Dropdown(choices=cols)

def dem_mods_change(value: str) -> Tuple[gr.Markdown, gr.DataFrame, gr.Textbox]:
    'Crop to Extent', 'Clip with Mask'
    if value == 'Crop to Extent':
        return gr.Markdown(visible=True), gr.DataFrame(visible=True), gr.Textbox(visible=False)
    elif value == 'Clip with Mask':
        return gr.Markdown(visible=False), gr.DataFrame(visible=False), gr.Textbox(visible=True)
    else:
        return gr.Markdown(visible=False), gr.DataFrame(visible=False), gr.Textbox(visible=False)

def write_mifn(dem, mifn, strm, spatial_units, flow_file, subtract_baseflow, rowcols_from_flowfile, flow_id, flow_params, flow_baseflow, vdt, is_database, num_iterations,
                                                    meta_file, convert_cfs_to_cms, x_distance, q_limit, lu_raster, is_lu_same_as_dem, mannings_table, direction_distance, slope_distance, low_spot_distance, low_spot_is_meters,
                                                    low_spot_use_box, box_size, find_flat, low_spot_find_flat_cutoff, degree_manip, degree_interval, Str_Limit_Val, UP_Str_Limit_Val, row_start, row_end, use_prev_d_4_xs,
                                                    weight_angles, man_n, adjust_flow, bathy_alpha, bathy_file, id_flow_file, omit_outliers, wse_search_dist, wse_threshold, wse_remove_three,
                                                    specify_depth, twd_factor, only_streams, use_ar_top_widths, flood_local, depth_map, flood_map, velocity_map, wse_map, fs_bathy_file,
                                                    da_flow_param,bathy_method,bathy_x_max_depth, bathy_y_shallow,fs_bathy_smooth_method, bathy_twd_factor) -> None:
    """
    Write the main input file
    """
    if not mifn:
        gr.Warning('Specify the main input file')
        return
    if not dem:
        gr.Warning('Specify the DEM')
        return
    
    # Format path strings
    mifn = _format_files(mifn)
    dem = _format_files(dem)
    
    if os.path.exists(mifn):
        gr.Warning(f'Overwriting main input file: {mifn}')
    
    if not isinstance(meta_file, bool) and os.path.exists(meta_file):
        gr.Warning(f'This will overwrite: {meta_file}')

    with open(mifn, 'w') as f:
        _warn_DNE('DEM', dem)
        _check_type('DEM', dem, ['.tif'])
        _write(f,'DEM_File',dem)

        f.write('\n')
        _write(f,'# AutoRoute Inputs')
        f.write('\n')

        if strm:
            strm = _format_files(strm)
            _warn_DNE('Stream File', strm)
            _check_type('Stream Fil', strm, ['.tif'])
            _write(f,'Stream_File',strm)

        _write(f,'Spatial_Units',spatial_units)

        if flow_file:
            flow_file = _format_files(flow_file)
            _warn_DNE('Flow File', flow_file)
            _check_type('Flow File', flow_file, ['.txt'])
            _write(f,'Flow_RAPIDFile',flow_file)

            if rowcols_from_flowfile:
                _write(f,'RowCol_From_RAPIDFile')
            if not flow_id:
                gr.Warning('Flow ID is not specified!!')
            else:
                _write(f,'RAPID_Flow_ID', flow_id)
            if not flow_params:
                gr.Warning('Flow Params are not specified!!')
            else:
                _write(f,'RAPID_Flow_Param', " ".join(flow_params))
            if subtract_baseflow:
                if not flow_baseflow:
                    gr.Warning('Base Flow Parameter is not specified, not subtracting baseflow')
                else:
                    _write(f,'RAPID_BaseFlow_Param',flow_baseflow)
                    _write(f,'RAPID_Subtract_BaseFlow')

        if vdt:
            vdt = _format_files(vdt)
            _check_type('VDT', vdt, ['.txt'])
            if is_database:
                    _write(f,'Print_VDT_Database',vdt)
                    _write(f,'Print_VDT_Database_NumIterations',num_iterations)
            else:
                _write(f,'Print_VDT',vdt)

        if meta_file:
            meta_file = _format_files(meta_file)
            _write(f,'Meta_File',meta_file)

        if convert_cfs_to_cms: _write(f,'CONVERT_Q_CFS_TO_CMS')

        _write(f,'X_Section_Dist',x_distance)
        _write(f,'Q_Limit',q_limit)
        _write(f,'Gen_Dir_Dist',direction_distance)
        _write(f,'Gen_Slope_Dist',slope_distance)
        _write(f,'Weight_Angles',weight_angles)
        _write(f,'Use_Prev_D_4_XS',use_prev_d_4_xs)
        _write(f,'ADJUST_FLOW_BY_FRACTION',adjust_flow)
        if Str_Limit_Val: _write(f,'Str_Limit_Val',Str_Limit_Val)
        if UP_Str_Limit_Val: _write(f,'UP_Str_Limit_Val',UP_Str_Limit_Val)
        if row_start: _write(f,'Layer_Row_Start',row_start)
        if row_end: _write(f,'Layer_Row_End',row_end)

        if degree_manip > 0 and degree_interval > 0:
            _write(f,'Degree_Manip',degree_manip)
            _write(f,'Degree_Interval',degree_interval)

        if lu_raster:
            lu_raster = _format_files(lu_raster)
            _warn_DNE('Land Use', lu_raster)
            _check_type('Land Use', lu_raster, ['.tif'])
            if is_lu_same_as_dem:
                _write(f,'LU_Raster_SameRes',lu_raster)
            else:
                _write(f,'LU_Raster',lu_raster)
            if not mannings_table:
                gr.Warning('No mannings table for the Land Use raster!')
            else:
                mannings_table = _format_files(mannings_table)
                _write(f,'LU_Manning_n',mannings_table)
        else:
            _write(f,'Man_n',man_n)

        if low_spot_distance:
            if low_spot_is_meters:
                _write(f,'Low_Spot_Dist_m',low_spot_distance)
            else:
                _write(f,'Low_Spot_Range',low_spot_distance)
            if low_spot_use_box:
                _write(f,'Low_Spot_Range_Box')
                _write(f,'Low_Spot_Range_Box_Size',box_size)
        
        if find_flat:
            if low_spot_find_flat_cutoff:
                _write(f,'Low_Spot_Find_Flat')
                _write(f,'Low_Spot_Range_FlowCutoff',low_spot_find_flat_cutoff)
            else:
                gr.Warning('Low Spot Range cutoff was not defined')

        if bathy_file:
            bathy_file = _format_files(bathy_file)
            _write(f,'BATHY_Out_File',bathy_file)
            _write(f,'Bathymetry_Alpha',bathy_alpha)

            if bathy_method == 'Parabolic':
                _write(f,'Bathymetry_Method',0)
            elif bathy_method == 'Left Bank Quadratic':
                _write(f,'Bathymetry_Method', 1)
                _write(f,'Bathymetry_XMaxDepth',bathy_x_max_depth)
                _write(f,'Bathymetry_YShallow',bathy_y_shallow)
            elif bathy_method == 'Right Bank Quadratic':
                _write(f,'Bathymetry_Method', 2)
                _write(f,'Bathymetry_XMaxDepth',bathy_x_max_depth)
                _write(f,'Bathymetry_YShallow',bathy_y_shallow)
            elif bathy_method == 'Double Quadratic':
                _write(f,'Bathymetry_Method', 3)
                _write(f,'Bathymetry_XMaxDepth',bathy_x_max_depth)
                _write(f,'Bathymetry_YShallow',bathy_y_shallow)
            elif bathy_method == 'Trapezoidal':
                _write(f,'Bathymetry_Method', 4)
                _write(f,'Bathymetry_XMaxDepth',bathy_x_max_depth)
            else: _write(f,'Bathymetry_Method', 5)

            if da_flow_param: _write(f, 'RAPID_DA_or_Flow_Param',da_flow_param)

        f.write('\n')
        _write(f,'# FloodSpreader Inputs')
        f.write('\n')

        if id_flow_file:
            id_flow_file = _format_files(id_flow_file)
            _warn_DNE('ID Flow File', id_flow_file)
            _check_type('ID Flow File', id_flow_file, ['.txt','.csv'])
            _write(f,'Comid_Flow_File',id_flow_file)

        if omit_outliers == 'Flood Bad Cells':
            _write(f,'Flood_BadCells')
        elif omit_outliers == 'Use AutoRoute Depths':
            _write(f,'FloodSpreader_Use_AR_Depths')
        elif omit_outliers == 'Smooth Water Surface Elevation':
            _write(f,'FloodSpreader_SmoothWSE_SearchDist',wse_search_dist)
            _write(f,'FloodSpreader_SmoothWSE_FractStDev',wse_threshold)
            _write(f,'FloodSpreader_SmoothWSE_RemoveHighThree',wse_remove_three)
        elif omit_outliers == 'Use AutoRoute Depths (StDev)':
            _write(f,'FloodSpreader_Use_AR_Depths_StDev')
        elif omit_outliers == 'Specify Depth':
            _write(f,'FloodSpreader_SpecifyDepth',specify_depth)

        _write(f,'TopWidthDistanceFactor',twd_factor)
        if only_streams: _write(f,'FloodSpreader_JustStrmDepths')
        if use_ar_top_widths: _write(f,'FloodSpreader_Use_AR_TopWidth')
        if flood_local: _write(f,'FloodLocalOnly')

        if depth_map:
            depth_map = _format_files(depth_map)
            _check_type('Depth Map',depth_map,['.tif'])
            _write(f,'OutDEP',depth_map)

        if flood_map:
            flood_map = _format_files(flood_map)
            _check_type('Flood Map',flood_map,['.tif'])
            _write(f,'OutFLD',flood_map)

        if velocity_map:
            velocity_map = _format_files(velocity_map)
            _check_type('Velocity Map',velocity_map,['.tif'])
            _write(f,'OutVEL',velocity_map)

        if wse_map:
            wse_map = _format_files(wse_map)
            _check_type('WSE Map',wse_map,['.tif'])
            _write(f,'OutWSE',wse_map)

        if fs_bathy_file and bathy_file: 
            fs_bathy_file = _format_files(fs_bathy_file)
            _check_type('FloodSpreader Generated Bathymetry',fs_bathy_file,['.tif'])
            _write(f,'FSOutBATHY', fs_bathy_file)
            if fs_bathy_smooth_method == 'Linear Interpolation':
                _write(f,'Bathy_LinearInterpolation')
            elif fs_bathy_smooth_method == 'Inverse-Distance Weighted':
                _write(f,'BathyTopWidthDistanceFactor', bathy_twd_factor)
       
async def run_exe(exe: str, mifn: str) -> None:
    """
    Run the executable, printing to terminal
    """
    if not exe or not mifn: return

    exe = _format_files(exe.strip())
    mifn = _format_files(mifn.strip())
    if not os.path.exists(exe):
        gr.Warning('AutoRoute Executable not found')
        return
    if not os.path.exists(mifn):
        gr.Warning('The Main Input File does not exist')
        return
    
    exe  = _prepare_exe(exe)

    process = await asyncio.create_subprocess_shell(f'echo "a" | {exe} {mifn}', # We echo a dummy input in so that AutoRoute can terminate if some input is wrong
                                                    stdout=asyncio.subprocess.PIPE,
                                                    stderr=asyncio.subprocess.PIPE)
    while True:
        line = await process.stdout.readline()
        if not line: 
            print('Program finished')
            break
        print(line.decode().strip())

    await process.communicate()
    
def _prepare_exe(exe: str) -> str:
    if SYSTEM == "Windows":
        return exe
    if SYSTEM == "Darwin" or SYSTEM == "Linux": # MAC
        slashes = exe.count('/')
        if slashes == 0:
            return './' + exe
        if slashes == 1 and exe[0] == '/':
            return '.' + exe
    return exe

def PrepareDEM(dem: str, 
               output: str, 
               mode: str, 
               extent: pd.DataFrame | None, 
               clip: str | None) -> None:
    """
    Prepare a DEM for use in AutoRoute / FloodSpreader programs. This function will:
        1) merge all DEMs into one file (if applicable)
        2) Clip DEM to an extent using the nearest cells
        3) Clip DEM to a region using the nearest cells

    Parameters
    ----------
    dem : any
        dem is either a path to a dem, or a directory containg .tif files. 
    extent: any, optional
        If specified, a list or tuple in the format of (minx, miny, maxx, maxy) is expected. 
        The DEM is clipped as close as possible to the extent
    clip : any, optional
        If specified, a file (.shp, .gpkg, .parquet) is expected. 
        The DEM is clipped as close as possible to the extent of the feature
    """
    # check the inputs
    if not dem or not output:
        gr.Warning("Please input both a DEM (or folder of DEMs) and the output file")
        return
    dem = _format_files(dem)
    output = _format_files(output)
    
    if mode == 'Crop to Extent':
        try:
            extent = [float(i) for i in extent.iloc[0].tolist()]
        except ValueError:
            gr.Warning("Please enter only numbers in each field for the Output Extent")
            return
        if extent[0] >= extent[2] or extent[1] >= extent[3]:
            gr.Warning(f"Invlaid extent {extent}")
            return
    elif mode == 'Clip with Mask':
        if clip: 
            clip = _format_files(clip)
        else:
            gr.Warning('Please specify a clipping file')
            return
        if not os.path.exists(clip) or not clip.endswith(('.shp','.gpkg','.parquet','.geoparquet')):
            gr.Warning(f"The provided clipping file '{clip}' does not exist or is unsupported")
            return
        if clip.endswith(('.shp','.gpkg')):
            extent = gpd.read_file(clip).total_bounds
        else:
            extent = gpd.read_parquet(clip).total_bounds
        
    if os.path.isfile(dem):
        dems = [dem]
    elif os.path.isdir(dem):
        dems = [f for f in glob.glob(os.path.join(dem, '*.tif')) if not os.path.basename(f).startswith('.')] # Exclude dot files that may be generated on some systems
        if len(dems) == 0:
            gr.Warning(f'No valid .tif files found in {dem}!')
            return
        if mode != 'None':
            new_dems = []
            for dem_file in dems:
                dataset = gdal.Open(dem_file)
                geotransform = dataset.GetGeoTransform()

                minx = geotransform[0]
                miny = geotransform[3] + geotransform[5] * dataset.RasterYSize
                maxx = geotransform[0] + geotransform[1] * dataset.RasterXSize
                maxy = geotransform[3]
                if (minx > extent[2] or maxx < extent[0] or miny > extent[3] or maxy < extent[1]):
                    new_dems.append(dem_file)
                dataset = None

            dems = new_dems
            if len(dems) == 0:
                gr.Warning(f'No .tif files were found in the specified extent in {dem}')
                return
    else:
        gr.Warning(f"Can't tell if {dem} is a file or a folder...")
        return

    # Create an in-memory VRT (Virtual Dataset) for merging the GeoTIFFs
    vrt_options = gdal.BuildVRTOptions(resampleAlg='bilinear')
    vrt_dataset = gdal.BuildVRT('', dems, options=vrt_options)
    projection = vrt_dataset.GetProjection()
    creation_options = ["COMPRESS=DEFLATE", "PREDICTOR=3","BIGTIFF=YES"] # Compression is compatible with old GDAL, big tiff should be yes to allow for some big tiffs

    if mode == 'Crop to Extent':
        print(projection, extent)
        warp_options = gdal.WarpOptions(
            format='GTiff',
            dstSRS=projection,
            dstNodata=-999,
            outputBounds=extent,
            outputType=gdal.GDT_Float32,  
            multithread=True, 
            creationOptions = creation_options
            )          
    elif mode == 'Clip with Mask':
        shp = ogr.Open(clip)
        layer = shp.GetLayer()
        name = layer.GetName()
        shp = None
        warp_options = gdal.WarpOptions(
            format='GTiff', 
            dstSRS=projection, 
            dstNodata=-999,
            cutlineDSName=clip, 
            cutlineLayer=name, 
            cropToCutline=True, 
            outputType=gdal.GDT_Float32, 
            multithread=True, 
            creationOptions = creation_options
            )
    else: # Just combine the tiles
        warp_options = gdal.WarpOptions(
            format='GTiff', 
            dstSRS=projection, 
            dstNodata=-999, 
            outputType=gdal.GDT_Float32, 
            multithread=True, 
            creationOptions = creation_options
            )

    try:
        print("Attempting warp...")
        gdal.Warp(output, vrt_dataset, options=warp_options)
    except RuntimeError as e:
        try:
            vrt_dataset = None
            disk, required = re.findall(r"(\d+)", str(e))
            gr.Warning(f"Need {_sizeof_fmt(int(required))}; {_sizeof_fmt(int(disk))} of space on this machine")
        except:
            gr.Warning(e)
        return
    
    # Clean up the VRT dataset
    vrt_dataset = None
    print("Finished 'Prepare DEM'")
    print()

def _sizeof_fmt(num:int) -> str:
    """
    Take in an int number of bytes, outputs a string that is human readable
    """
    for unit in ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB"):
        if abs(num) < 1024.0:
            return f"{num:3.1f} {unit}"
        num /= 1024.0
    return f"{num:.1f} YB"