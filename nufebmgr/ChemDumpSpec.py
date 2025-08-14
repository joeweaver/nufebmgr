class ChemDumpSpec:
    """
    Specification for chemical-related output variables in simulation dumps.

    This class centralizes the configuration of which chemical-related variables
    should be included in output files. Currently, only HDF5 dumps are supported;
    support for VTK dumps is planned.

    Parameters
    ----------
    outspec_h5 : str, optional
        Specification of which variables to include in HDF5 output. Must be one of:
        - ``"all"``: Include all available variables (default).
        - ``"conc"``: Include only ID, type, and spatial coordinates.
        - ``"reac"``: Include radius alongside ID, type, and spatial coordinates.
        - ``"custom"``: Include only the variables listed in ``custom_h5``.
    custom_h5 : list of str, optional
        List of variables to include when ``outspec_h5`` is set to ``"custom"``.

    Notes
    -----
    Rationale for defaults:
        The default output is ``"all"`` because it is easier to notice an excessively
        large dump file and trim output in subsequent runs than to re-run many
        simulations because essential data was missing.

    Raises
    ------
    ValueError
        If ``outspec_h5`` is not one of the valid options.

    Examples
    --------
    Include all variables in HDF5 output:

        >>> spec = ChemDumpSpec()
        >>> spec.hdf5_vars()
        ['conc', 'reac']

    Include only location variables:

        >>> spec = ChemDumpSpec(outspec_h5="conc")
        >>> spec.hdf5_vars()
        ['conc']

    Include custom variables (note that there is no checking on custom variable suitability:

        >>> spec = ChemDumpSpec(outspec_h5="custom", custom_h5=["banana", "reac"])
        >>> spec.hdf5_vars()
        ['banana', 'reac']
    """

    _VALID_H5 = ["all", "conc", "reac", "custom"]

    def __init__(self, outspec_h5: str = "all", custom_h5: list[str] | None = None):
        if outspec_h5 not in ChemDumpSpec._VALID_H5:
            raise ValueError(
                f'Unknown output specification "{outspec_h5}". Valid values are: {", ".join(ChemDumpSpec._VALID_H5)}')
        self.outspec_h5 = outspec_h5
        if custom_h5 is None:
            custom_h5 = []
        self.custom_h5 = custom_h5

    def hdf5_vars(self) -> list[str]:
        """
        Return the list of bug-related HDF5 variable names to dump.

        The exact list depends on the `outspec_h5` value provided during
        initialisation

        Returns
        -------
        list[str]
            Variable names to include in the HDF5 output.
        """

        if self.outspec_h5 == "all":
            return ['conc', 'reac']
        if self.outspec_h5 == "conc":
            return ['conc']
        if self.outspec_h5 == "reac":
            return ['reac']
        if self.outspec_h5 == "custom":
            return self.custom_h5
        else:
            raise ValueError(f'Unknown output specification {self.outspec_h5}. Valid values are: {", ".join(ChemDumpSpec._VALID_H5)}')