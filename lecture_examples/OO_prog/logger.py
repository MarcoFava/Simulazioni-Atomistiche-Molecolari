verbosity = "INFO"

def info(message):
    if verbosity == "INFO":
        print('INFO:', message)
    if verbosity == "WARNING":
        print('WARNING: ', message)