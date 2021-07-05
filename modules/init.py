import io

def readConfig(config_file):
    config={}
    with open(config_file,'r') as CONFIG:
        for line in CONFIG:
            line = line.strip()
            if line.startswith("#"):
                continue
            contents = line.split("=")
            if len(contents)==2:
                key = contents[0].strip()
                value = contents[1].strip().strip('"')
                config[key] = value
    return config