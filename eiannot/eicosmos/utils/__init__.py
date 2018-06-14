from eicore.external_process.snakemake_helper import loadPreCmd


def loader(programs=()):
    def in_loader(func):
        """This wrapper has to be used as decorator.
        It will create, if necessary, the necessary string for preloading the environment on the command line."""
        def function_wrapper(*args, **kwargs):
            cc = ""
            for program in programs:
                pre = kwargs["config"].get("load", dict()).get(program, "")
                if pre:
                    cc += "{} &&".format(pre)
            if cc:
                cc = "set +u && {}".format(cc)
                kwargs["load"] = cc

            return func(*args, **kwargs)
        return function_wrapper
    return in_loader


def define_output(func):

    def function_wrapper(*args, **kwargs):
        if "in_file" in kwargs and "suffix" in kwargs and "out_file" not in kwargs:
            kwargs["out_file"] = kwargs["out_file"] + kwargs["suffix"]
        return func(*args, **kwargs)
    return function_wrapper

