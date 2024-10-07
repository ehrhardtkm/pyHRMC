# validator architecture taken from Simmate
# Simmate is licensed undr a BSD license, which can be found in LICENSE.md


class Validator:
    @classmethod
    @property
    def name(cls):
        """
        A nice string name for the validator. By default it just returns the name
        of this class by default.
        """
        return cls.__name__
