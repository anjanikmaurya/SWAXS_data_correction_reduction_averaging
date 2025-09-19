import numpy as np
from pathlib import Path
import logging
from typing import List
logger = logging.getLogger('swaxs_pipeline')

def read_raw_detector_15_image(raw_file_path: Path, shape: List[int], log_pending: bool) -> np.ndarray:
    """
    Read raw detector file and return detector image array.
    
    Parameters
    ----------
    raw_file_path Path to the .raw detector file
    detector_type : Type of detector ('SAXS' or 'WAXS')
        
    Returns
    -------
    np.ndarray
        Detector image array with proper shape
    """
    # Log function call for first runs only
    if log_pending:
        logger.info(f"Called read_raw_detector_15_image()")
    data = np.fromfile(str(raw_file_path), dtype=np.int32).reshape(shape)

    
    if not raw_file_path.exists():
        raise FileNotFoundError(f"Raw detector file not found: {raw_file_path}")
    
    data = np.fromfile(str(raw_file_path), dtype=np.int32).reshape(shape)
    
    # Log data information for first runs only
    if log_pending:
        logger.info(f"  Raw data loaded: shape {data.shape}, dtype {data.dtype}")
        logger.info(f"  Data range: min={data.min()}, max={data.max()}, mean={data.mean():.2f}")
    
    return data