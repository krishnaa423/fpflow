#region modules
import logging 
#endregion

#region variables
log_mode = logging.INFO
logging.basicConfig(level=log_mode)
logger = logging.getLogger('AppLogger')
#endregion

#region functions
def get_logger():
    return logger
#endregion

#region classes
#endregion