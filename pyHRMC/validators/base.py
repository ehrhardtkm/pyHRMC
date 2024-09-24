#validator architecture taken from Simmate

class Validator:
    @classmethod
    @property
    def name(cls):
        """
        A nice string name for the validator. By default it just returns the name
        of this class by default.
        """
        return cls.__name__