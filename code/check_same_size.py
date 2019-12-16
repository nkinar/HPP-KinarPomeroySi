from get_size import length

def check_same_size(a, b, err=''):
    if length(a) != length(b):
        raise ValueError(err)
# DONE
