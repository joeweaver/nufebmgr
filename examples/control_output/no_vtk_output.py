from nufebmgr import NufebProject


def configure_project():
    with NufebProject() as prj:
        prj.use_seed(1701)
        prj.set_box(x=50,
                    y=50,
                    z=50)
        prj.add_taxon_by_template(name="basic_het",
                                  template="basic_heterotroph")
        prj.add_taxon_by_template(name="slow_het",
                                  template="slow_heterotroph")
        prj.add_taxon_by_template(name="small_het",
                                  template="small_heterotroph")

        prj.layout_uniform(nbugs=15)

        prj.set_composition({'basic_het':'33.33',
                             'slow_het':'33.33',
                             'small_het':'33.33',})

        prj.distribute_spatially_even()

        prj.set_track_abs()
        prj.enable_thermo_output(timestep=1)
        prj.disable_vtk_output()

        atom_in, inputscript = prj.generate_case()
        return atom_in, inputscript


if __name__ == '__main__':
    config, inputscript = configure_project()

    with open("atom.in", "w") as outf:
        outf.write(config)

    with open("inputscript.nufeb", "w") as outf:
        outf.write(inputscript)
