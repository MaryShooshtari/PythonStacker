import argparse
import os
import sys
import json
import glob
import src.jobSubmission.condorTools as ct

import src.arguments as arguments


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)

    arguments.add_settingfiles(parser)
    arguments.add_tmp_storage(parser)
    arguments.add_reweighting_option(parser)
    args, args_unknown = parser.parse_known_args()

    if len(args_unknown) > 0 or len(sys.argv) == 1:
        parser.print_help()
        parser.exit()
    return args


if __name__ == "__main__":
    args = parse_arguments()

    # loop processes and check which one has EFT, then submit it
    with open(args.processfile, 'r') as f:
        content = json.load(f)
        processes = content["Processes"]
        basedir = content["Basedir"]

    commandset = []
    for year in args.years:
        if year == "2016":
            continue
        for pname, process in processes.items():
            # if has EFT
            files = []
            for filebase in process["fileglobs"]:
                fileglob = os.path.join(basedir, filebase)
                # fileglob += f"*{args.years[0]}"
                fileglob += f"*{year}*base.root"
                #print(fileglob)
                true_files = glob.glob(fileglob)
                files.extend(true_files)

            cmds = []
            for filename in files:
                if process.get("hasBSM", 0) == 0:
                    continue
                for eventclass in range(14):
                    cmd = "python3 createBSMWeights.py"
                    cmd += f" -f {filename}"
                    cmd += f" -e {eventclass}"
                    cmd += f" -p {pname}"
                    cmd += f" --storage {args.storage}"
                    if process.get('doReweight', 0) != 0:
                        cmd += f" --pseudo {args.pseudo}"
                    cmds.append(cmd)

            if len(cmds) != 0:
                # cmds.append("echo 'nothing to do'")
                commandset.append(cmds)

    ct.submitCommandsetsAsCondorCluster("CreateBSMWeights", commandset, scriptfolder="Scripts/condor/")
