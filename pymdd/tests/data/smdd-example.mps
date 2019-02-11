NAME           smdd-example
ROWS
 N  obj     
 G  c1      
 G  c2      
 L  c3      
COLUMNS
    MARK0000  'MARKER'                 'INTORG'
    x1        obj                  4   c1                   1
    x1        c3                   1
    x2        obj                  3   c2                   1
    x2        c3                   1
    x3        obj                  2   c1                   1
    x3        c2                   1   c3                   1
    MARK0001  'MARKER'                 'INTEND'
RHS
    rhs       c1                   1   c2                   1
    rhs       c3                   2
BOUNDS
 UP bnd       x1                   1
 UP bnd       x2                   1
 UP bnd       x3                   1
ENDATA
