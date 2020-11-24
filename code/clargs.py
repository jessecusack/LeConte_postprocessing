import os
from argparse import ArgumentParser
from munch import munchify


def gen_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "-ld",
        "--leconte",
        dest="leconte",
        default=os.path.expanduser("~/Dropbox/LeConte"),
        help="path to LeConte Dropbox directory",
    )
    parser.add_argument(
        "-ad",
        "--adcp",
        dest="adcp",
        default=os.path.expanduser("~/Dropbox/LeConte_ADCP_final"),
        help="path to LeConte ADCP directory",
    )
    parser.add_argument(
        "-sd",
        "--save",
        dest="save",
        default="../proc",
        help="path to save processed data",
    )
    # parser.add_argument("-q", "--quiet",
    #                     action="store_false", dest="verbose", default=False,
    #                     help="don't print status messages to stdout")

    return parser


def check_args(args):
    # Check specified directories exist.

    if not os.path.exists(args.leconte):
        raise ValueError(
            "LeConte data directory does not exist: '{}'".format(args.leconte)
        )
    else:
        print("LeConte data path '{}' exists.".format(args.leconte))

    if not os.path.exists(args.adcp):
        raise ValueError("ADCP directory does not exist: '{}'".format(args.adcp))
    else:
        print("ADCP path '{}' exists.".format(args.adcp))

    if not os.path.exists(args.save):
        raise ValueError("Save directory does not exist: '{}'".format(args.save))
    else:
        print("Save path '{}' exists.".format(args.save))

        
def parse_check_args():
    # Combine parsing and checking.
    parser = gen_parser()
    args = munchify(vars(parser.parse_args()))
    check_args(args)
    return args