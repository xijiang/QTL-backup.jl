"""
    function readped(pf)
Read `ID`, `Sire`, and `Dam` from file `pf` into a DataFrame.
The DataFrame will have 4 columns:
- ID: recoded names in integer start from 0 for unknown
- Sire and dam columns: integer, 0 for unknown
- Name: the original name as `String`
"""
function readped(pf; skip = 0)
end

