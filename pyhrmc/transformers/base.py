#architecture from Simmate

class Transformation:
    @classmethod
    @property
    def name(cls):
        """
        A nice string name for the transformation. By default it just returns
        the name of this class.
        """
        return cls.__name__