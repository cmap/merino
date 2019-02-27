import yaml

def write_args_to_file(args, outfile):
    with open(outfile, "w") as f:
        print vars(args)
        yaml.dump(vars(args), f)