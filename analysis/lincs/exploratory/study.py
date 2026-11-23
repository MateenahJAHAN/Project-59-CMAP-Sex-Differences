"""Open the released model configuration against this run's verified matrices."""
from pathlib import Path
import sys

ENGINE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ENGINE / 'src/G'))
import gconfig
import ginputs


def open_study(work):
    work = Path(work).resolve()
    config = gconfig.load(str(ENGINE / 'config/G_defaults.json'))
    paths = {'inputs':str(work / 'model_inputs'), 'tensors_D':str(work / 'data/D/tensors'),
             'tensors_E':str(work / 'data/E/tensors'), 'persig_E':str(work / 'data/E/tensors'),
             'siginfo':str(work / 'sources/siginfo_beta.txt'),
             'prism':str(work / 'sources/prism_secondary_dose_response_curve_parameters.csv')}
    return config, ginputs.Inputs(str(ENGINE), config, paths, strict=True)
