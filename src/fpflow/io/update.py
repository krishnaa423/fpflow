#region modules
#endregion

#region variables
#endregion

#region functions
def update_dict(dest: dict, source: dict) -> dict:
    if source is None or source=={}:
        return dest 

    for source_key, source_value in source.items():
        if source_key in dest.keys() and isinstance(source_value, dict):
            dest[source_key].update(source_value)
        else:
            dest[source_key] = source_value

    return dest
#endregion

#region classes
#endregion