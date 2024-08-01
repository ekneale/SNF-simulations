import numpy as np
import sample

def detector_signal_calc(total_spec,efficiency):
    events = np.array(sample.sample(total_spec, 1000))

    detectable = events[np.where(events > 1800)]

    detector_sig = efficiency * len(detectable)

    return detector_sig



