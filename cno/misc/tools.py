__all__ = ["CNOError"]
    

class CNOError(Exception):
    """A simple exception related to CellNOpt"""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

