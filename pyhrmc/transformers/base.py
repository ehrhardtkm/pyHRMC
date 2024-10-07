# architecture from Simmate
# Simmate is licensed under a BSD license, which can be found in LICENSE.md


class Transformation:
    @classmethod
    @property
    def name(cls):
        """
        A nice string name for the transformation. By default it just returns
        the name of this class.
        """
        return cls.__name__
