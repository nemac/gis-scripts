import os
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.windows import Window
import numpy as np

# === PERFORMANCE SETTINGS ===
os.environ["GDAL_NUM_THREADS"] = "ALL_CPUS"
os.environ["GDAL_CACHEMAX"] = "8192"  # Increased to 8GB

# === USER SETTINGS ===
input_raster = r"D:\AWS_download\Mosaicked\BH4_BH5_mosaic.tif"
output_raster = r"D:\AWS_download\Mosaicked\NA_4m_proj_30m.tif"
target_crs = "ESRI:102003"
target_res = 30

def reproject_optimized(in_raster, out_raster, target_crs, target_res):
    with rasterio.Env(GDAL_NUM_THREADS="ALL_CPUS", GDAL_CACHEMAX=8192):
        with rasterio.open(in_raster) as src:
            print(f"📘 Source: {src.width}x{src.height}, {src.count} bands")
            print(f"🧭 Target CRS: {target_crs}, Resolution: {target_res}m")
            
            # Calculate transform
            transform, width, height = calculate_default_transform(
                src.crs, target_crs, src.width, src.height, *src.bounds,
                resolution=target_res
            )
            
            # Update metadata with compression and tiling
            kwargs = src.meta.copy()
            kwargs.update({
                'crs': target_crs,
                'transform': transform,
                'width': width,
                'height': height,
                'compress': 'lzw',      # Add compression
                'tiled': True,           # Enable tiling
                'blockxsize': 512,       # Tile size
                'blockysize': 512,
                'BIGTIFF': 'YES'         # Support large files
            })
            
            # Create output file
            with rasterio.open(out_raster, 'w', **kwargs) as dst:
                # Reproject ALL bands at once (key optimization!)
                for band_idx in range(1, src.count + 1):
                    print(f"Processing band {band_idx}/{src.count}")
                    
                    reproject(
                        source=rasterio.band(src, band_idx),
                        destination=rasterio.band(dst, band_idx),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=transform,
                        dst_crs=target_crs,
                        resampling=Resampling.bilinear,
                        num_threads=os.cpu_count()
                    )
            
            print(f"✅ Complete: {out_raster}")

# === RUN ===
reproject_optimized(input_raster, output_raster, target_crs, target_res)
