def rate_of_spread(sa_vol_ratio, fuel_load, bed_depth, dead_extinction_moisture, moisture_content, wind, tan_slope,
                    heat_content = 8000, mineral_content = .0555,
                    effective_mineral_content = .010, particle_density = 32):
    """ Basic Rothermel Model taken from:
    https://www.fs.fed.us/rm/pubs_series/rmrs/gtr/rmrs_gtr371.pdf
    :param sa_vol_ratio: surface to volume ratio in ft^-1
    :param fuel_load: fuel load in lb/ft^2
    :param bed_depth: fuel bed depth in ft
    :param dead_extinction_moisture: dead fuel moisture of extinction (fraction)
    :param moisture_content: moisture content (fraction)
    :param wind: wind velocity at midflame height in ft/min)
    :param tan_slope: slope steepness, maximum (fraction)
    :param heat_content: low heat content, 8000 Btu/lb by default
    :param mineral_content: total mineral content (fraction) 0.0555 lb mineral /lb wood by default
    :param effective_mineral_content: effective mineral content (fraction) .010 by default
    :param particle_density: Oven-dry particle density 32 lb/ft^3 by default
    :return: rate of spread of the fire in ft/min
    """
    preignition = 250 + 1116 * moisture_content

    heating_number = np.exp(-138 / sa_vol_ratio)

    bulk_density = fuel_load / bed_depth
    packing_ratio = bulk_density / particle_density
    optimal_packing = 3.348 * sa_vol_ratio ** -.8189

    slope_factor = 5.275 * packing_ratio ** -.3 * tan_slope**2

    C = 7.74 * np.exp(-.133 * sa_vol_ratio**.55)
    B = .02526 * sa_vol_ratio ** .54
    E = .715 * np.exp(sa_vol_ratio * -3.59e-4)
    wind_factor = C * wind**B * (packing_ratio / optimal_packing) ** E

    propegating_flux = (192 +.2595 * sa_vol_ratio) **-1 * np.exp((.792 + .681 * sa_vol_ratio ** .5) * (packing_ratio + .1))

    mineral_dampening = min(1.0, .174 * effective_mineral_content**-.19)

    rm = min(1.0, moisture_content / dead_extinction_moisture)
    moisture_dampening = 1 - 2.59 * rm + 5.11 * rm**2 - 3.52 * rm**3

    net_fuel_load = fuel_load * (1 - mineral_content)

    max_reaction = sa_vol_ratio ** 1.5 * (495 + .0594 * sa_vol_ratio**1.5)**-1

    A = 133 * sa_vol_ratio ** -.7913

    optimal_reaction = max_reaction * (packing_ratio / optimal_packing) ** A * np.exp(A * (1 - packing_ratio / optimal_packing))

    reaction_intensity = optimal_reaction * net_fuel_load * heat_content * moisture_dampening * mineral_dampening

    rate_of_spread = reaction_intensity * propegating_flux * (1 + wind_factor + slope_factor) / (bulk_density * heating_number * preignition)

    return rate_of_spread
