complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-':'-'}
def rc(seq):
    return "".join(complement.get(base, base) for base in reversed(seq))