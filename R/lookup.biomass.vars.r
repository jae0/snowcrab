
  lookup.biomass.vars = function() {

    # vars in ("name", "filter") pairs:
      vars= matrix( c(
            "totmass.male",             "male",
            "totmass.male.com",         "m.com",
            "totmass.male.ncom",        "m.ncom",
            "totmass.female",           "female",
            "totmass.female.berried",   "f.berried",
            "totmass.female.primiparous","primiparous",
            "totmass.female.multiparous", "multiparous",
            "totmass.female.mat",       "f.mat",
            "totmass.female.imm",       "f.imm",
            "totmass.male.mat",         "m.mat",
            "totmass.male.imm",         "m.imm",
            "totmass.female.soft",      "f.soft",
            "totmass.male.hard",      "m.hard",
            "totmass.male.soft",      "m.soft",
            "totmass.female.hard",      "f.hard",
            "totmass.male.com.CC1to2",      "m.CC1to2",
            "totmass.male.com.CC3to4",      "m.CC3to4",
            "totmass.male.com.CC1",         "m.CC1",
            "totmass.male.com.CC2",         "m.CC2",
            "totmass.male.com.CC3",         "m.CC3",
            "totmass.male.com.CC4",         "m.CC4",
            "totmass.male.com.CC5",         "m.CC5",
            "totmass.female.CC1to2",      "f.CC1to2",
            "totmass.female.CC3to4",      "f.CC3to4",
            "totmass.female.CC3",      "f.CC3",
            "totmass.female.CC4",      "f.CC4",
            "totmass.female.CC5",         "f.CC5",
            "mi123.mass",                  "mi123",
            "mi4.mass",                    "mi4",
            "mi5.mass",                    "mi5",
            "mi6.mass",                    "mi6",
            "mi7.mass",                    "mi7",
            "mi8.mass",                    "mi8",
            "mi9.mass",                    "mi9",
            "mi10.mass",                    "mi10",
            "mi11.mass",                    "mi11",
            "mi12.mass",                    "mi12",
            "fi1234.mass",                  "fi1234",
            "fi5.mass",                    "fi5",
            "fi6.mass",                    "fi6",
            "fi7.mass",                    "fi7",
            "fi8.mass",                    "fi8",
            "fi9.mass",                    "fi9",
            "fi10.mass",                   "fi10",
            "R0.mass",                  "R0",
            "R1.mass",                  "R1",
            "R2.mass",                  "R2",
            "R3.mass",                  "R3",
            "R4.mass",                  "R4",
            "R5p.mass",                 "R5p",
            "male.small.mass",          "male.small",
            "male.large.mass",          "male.large",
            "female.small.mass",          "female.small",
            "female.large.mass",          "female.large",
             "dwarf.mass",               "m.dwarf",
            "pre.recruit.mass",         "pre.recruit",
            "totmass.male.skip.moulter",        "skip.moulter"
            ), ncol=2, byrow=T)

    return ( vars )    

  }


