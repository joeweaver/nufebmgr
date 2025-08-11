from nufebmgr import NufebProject


def configure_project():
    with NufebProject() as prj:
        prj.use_seed(1701)
        prj.set_box(x=50,
                    y=50,
                    z=50)
        prj.add_taxon_by_template(name="attacker",
                                  template="basic_heterotroph")
        prj.add_taxon_by_template(name="vuln",
                                  template="slow_heterotroph")
        prj.add_taxon_by_template(name="immune",
                                  template="small_heterotroph")

        prj.layout_poisson(5)

        prj.set_composition({'attacker_slow':'33.33',
                             'vuln':'33.33',
                             'immune':'33.33',})

        prj.distribute_spatially_even()

        prj.set_track_abs()
        prj.enable_thermo_output(timestep=1)
        prj.disable_hdf5_output()
        prj.add_lysis_group_by_json('vuln_intoxicated',{'name':'vuln_intoxicated','releases':'sub','rate':'2e-3','percent':'0.2'})

        prj.arm_t6ss(taxon="attacker",effector="toxin_a",harpoon_len=1.3e-6,cooldown=100)
        prj.vuln_t6ss(taxon="immune",effector="toxin_a",prob=1,to_group="vuln_intoxicated")
        prj.limit_biofilm_height(2)
        prj.set_runtime(4*7*24*60*60)

        atom_in, inputscript = prj.generate_case()
        return atom_in, inputscript


if __name__ == '__main__':
    config, inputscript = configure_project()

    with open("atom.in", "w") as outf:
        outf.write(config)

    with open("inputscript.nufeb", "w") as outf:
        outf.write(inputscript)
