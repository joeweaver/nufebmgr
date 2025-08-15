from nufebmgr import NufebProject


def configure_project():
    with NufebProject() as prj:
        prj.use_seed(1701)
        prj.set_box(x=120,
                    y=120,
                    z=50)
        prj.enable_inferring_substrates()
        prj.add_taxon_by_template(name="basic_het",
                                  template="basic_heterotroph")
        prj.add_taxon_by_template(name="slow_het",
                                  template="slow_heterotroph")
        prj.add_taxon_by_template(name="small_het",
                                  template="small_heterotroph")
        prj.set_taxa_groups({'basic_het': '4', 'slow_het': '6', 'small_het': '7'})

        prj.layout_poisson(radius=8)


        prj.distribute_even_strips("horizontal", noise=20)

        prj.set_track_abs()
        prj.enable_thermo_output(timestep=1)
        atom_in, inputscript = prj.generate_case()
        return atom_in, inputscript


if __name__ == '__main__':
    config, inputscript = configure_project()

    with open("atom.in", "w") as outf:
        outf.write(config)

    with open("inputscript.nufeb", "w") as outf:
        outf.write(inputscript)
