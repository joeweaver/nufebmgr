from nufebmgr import NufebProject

def configure_project():
    with NufebProject() as prj:
        prj.use_seed(42)
        prj.set_box(x=50,
                    y=50,
                    z=50)
        prj.enable_inferring_substrates()
        prj.add_taxon_by_jsonfile("example.json")

        prj.layout_uniform(nbugs=20)

        prj.set_composition({'basic_heterotroph':'25',
                             'slow_heterotroph':'50',
                             'small_heterotroph':'25'})
        prj.set_taxa_groups({'basic_het': '4', 'slow_het': '6', 'small_het': '7'})
        prj.distribute_spatially_even()

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


