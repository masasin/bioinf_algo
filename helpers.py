#!/usr/bin/env python
def read_file(path_in, line_types='s'):
    '''
    `line_types` is a character sequence, with `s` for string, `i` for int, and
    `f` for float.

    >>> genome, number = read_file(some_path, 'si')  # doctest: +SKIP

    '''
    with open(path_in, 'r') as fin:
        for line_type in line_types:
            line = fin.readline().strip()
            if line_type == 's':
                yield line
            elif line_type == 'i':
                yield int(line)
            elif line_type == 'f':
                yield float(line)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
