from snf_simulations import scale, load_spec

# combining load and scale functions


def load_equal_scaled(
    data, max_E, named, isotope, m, mr, half_life_yrs, removal_time, min_E=0
):
    spec = load_spec.load_equal(
        named, isotope, data[:, 0], data[:, 1], data[:, 2], max_E, min_E
    )
    spec_scaled = scale.scale(spec, m, mr, half_life_yrs, removal_time)
    spec_scaled.SetTitle(isotope)
    return spec_scaled
