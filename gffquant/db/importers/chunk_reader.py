""" chunked file reader """


def get_lines(stream, bufsize=4000000):
    tail = ""
    while 1:
        chunk = "".join((tail, stream.read(bufsize).decode()))
        if not chunk:
            break
        chunk = chunk.split("\n")
        if len(chunk) > 1:
            chunk, tail = chunk[:-1], chunk[-1]
        for line in chunk:
            if line[0] != "#":
                yield line
