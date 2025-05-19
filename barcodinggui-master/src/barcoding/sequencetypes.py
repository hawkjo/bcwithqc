from enum import IntEnum


class SequencingMode(IntEnum):
    singleEnd = 1
    pairedEnd = 2


class SequenceType(IntEnum):
    def __new__(cls, value, name):
        obj = int.__new__(cls, value)
        obj._value_ = value
        return obj

    def __init__(self, value, name):
        self.display_name = name

    randomBarcode = 0, "random sequence"
    constantRegion = 1, "constant region"
    barcodeList = 2, "barcode whitelist"
