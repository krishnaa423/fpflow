#region modules
import functools
import os 
#endregion

#region variables
#endregion

#region functions
def change_dir(method):
    """Decorator for class methods: 
    temporarily change to self.dest_dir during execution."""
    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        os.chdir(self.dest_dir)
        result =  method(self, *args, **kwargs)
        os.chdir(self.current_dir)

        return result 
    return wrapper
#endregion

#region classes
#endregion