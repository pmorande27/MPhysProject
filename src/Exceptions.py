class ChiralExceptions(Exception):
    def __init__(self,message) -> None:
        self.message = message
        super().__init__(message)
class ThermalizationException(ChiralExceptions):
    def __init__(self,message) -> None:
        self.message = message
        super().__init__(message)

class CalibrationException(ChiralExceptions):
    def __init__(self,message) -> None:
        self.message = message
        super().__init__(message)
