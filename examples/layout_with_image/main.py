from nufebmgr import NufebProject


def configure_project():
    with NufebProject() as prj:
        prj.use_seed(1701)
        prj.set_box(x=100,
                    y=100,
                    z=50)
        prj.enable_inferring_substrates()
        prj.add_taxon_by_jsonfile("example.json")

        mappings = {"FF1B9E77": 'basic_heterotroph',
                    "FFD95F02": 'slow_heterotroph',
                    "FF7570B3": 'small_heterotroph'}

        prj.simple_image_layout('layout.png', mappings)
        prj.set_track_abs()
        prj.enable_thermo_output(timestep=1)
        prj.run_for_N_steps(int(6*60*60/900)) # 6 hours assuming 900s biological timestep

        atom_in, inputscript = prj.generate_case()
        return atom_in, inputscript


if __name__ == '__main__':
    config, inputscript = configure_project()

    with open("atom.in", "w") as outf:
        outf.write(config)

    with open("inputscript.nufeb", "w") as outf:
        outf.write(inputscript)
