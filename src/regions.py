"""Study regions as lists of county/parish FIPS codes (Census cartographic boundaries).

`gulf` = counties/parishes that front the Gulf of Mexico or its bays/sounds (TX-FL).
This list was assembled by hand from geography, not taken from an official NOAA file.
Inland counties that NOAA counts as coastal-watershed counties are NOT included. Review before citing.
"""
REGIONS = {
    "houston": {"TX": ["48201"]},
    "gulf": {
        "TX": "48061 48489 48261 48273 48355 48409 48007 48391 48057 48321 48039 48167 48071 48245 48361 48201 48239".split(),
        "LA": "22023 22113 22045 22101 22109 22057 22051 22075 22087 22071 22103 22019".split(),
        "MS": "28045 28047 28059".split(),
        "AL": "01097 01003".split(),
        "FL": "12033 12113 12091 12131 12005 12045 12037 12129 12065 12123 12029 12075 12017 12053 12101 12103 12057 12081 12115 12015 12071 12021 12087".split(),
    },
}
