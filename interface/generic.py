def split_geom_and_grid(geom):
    # split data into two sparate lists
    # as one will process many log files and want only 1 geometry but the full
    # ims grid
    g = geom
    ims_grid = []
    geom = []
    for el in g:
        if "Bq" in el['label']:  # it is a bq atom -> ims grid
            ims_grid.append(el)
        else:
            geom.append(el)
    return geom, ims_grid