import json

from src.configuration.Uncertainty import Uncertainty


def load_uncertainties(jsonfile, typefilter=None, namefilter=None, allowflat=True) -> dict:
    with open(jsonfile, 'r') as f:
        all_data = json.load(f)
    ret = dict()
    for key, val in all_data.items():
        if typefilter and val.get("type", "flat") != typefilter:
            continue
        if namefilter and val["name"] != namefilter:
            continue
        if not allowflat and val.get("type", "flat") == "flat":
            continue

        ret[key] = Uncertainty(key, val)
    return ret


def load_channels(channelfile) -> dict:
    ret = dict()
    with open(channelfile, 'r') as f:
        data = json.load(f)
    for key, val in data.items():
        channel_tmp = Channel(val, data)
        if channel_tmp.isSubchannel:
            continue
        ret[key] = channel_tmp
    return ret


# TODO: fix init
class Channel:
    def __init__(self, channelinfo, full_channelfile) -> None:
        self._selection = channelinfo["selection"]
        self._isSubchannel = bool(channelinfo.get("isSubchannel", 0))

        # subchannel structure
        subchannels = channelinfo.get("subchannels", [])
        self.subchannels = {subchannel: Channel(full_channelfile[subchannel], full_channelfile) for subchannel in subchannels}

        # build process ignore list
        self.ign_processes = channelinfo.get("ignore_processes")

        # load other config options

    def is_process_excluded(self, process: str):
        return process in self.ign_processes

    def produce_masks(self, tree):
        # use produce aliases to produce boolean masks for the subchannels
        aliases, keys = self.produce_aliases()
        masks = tree.arrays(keys, cut=self.selection, aliases=aliases)
        return masks, keys

    def produce_aliases(self) -> tuple[dict, list]:
        # for the chosen channel, load the info
        aliases: dict = {}
        for name, info in self.subchannels.items():
            aliases[name] = info.selection

        alias_names = list(aliases.keys())
        return aliases, alias_names

    @property
    def selection(self) -> str:
        return self._selection

    @selection.setter
    def selection(self, value):
        if not isinstance(value, str):
            raise ValueError("selection must be a str value")
        self._selection = value

    @property
    def isSubchannel(self) -> bool:
        return self._isSubchannel

    @isSubchannel.setter
    def isSubchannel(self, value):
        if not isinstance(value, bool):
            raise ValueError("isSubchannel must be a bool value")
        self._isSubchannel = value


# might replace with a class just keeping track of a few arrays.
# TODO: start from json file or at least dict
class Process:
    def __init__(self) -> None:
        self._isSignal = False
        self.name = ""

    @property
    def isSignal(self):
        return self._isSignal

    @isSignal.setter
    def isSignal(self, value):
        if not isinstance(value, bool):
            raise ValueError("isSignal must be a boolean value")
        self._isSignal = value