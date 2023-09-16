use ExtUtils::MakeMaker;
     WriteMakefile(
             'NAME'    => 'p5namot',                 # Name of package
             'INC' => '-I/home/tpascal/codes/namot2/source/namot-2.2.0-pre4/src/ -I/home/tpascal/codes/namot2/source/namot-2.2.0-pre4/src//../   -Wformat=0 -Wno-format-security -Wno-int-to-pointer-cast',  
             'LIBS'    => ['-lutil -ldl -lpthread -lpng -lm -lXm -lXt -lXext -lgsl -lgslcblas -lX11  -lSM -lICE'],  # Name of custom libraries
             'OBJECT'  => 'p5namot_wrap.o p5namot-namot_wrap.o p5namot-Backbone.o p5namot-dat_man.o p5namot-main.o p5namot-read_amber.o p5namot-rms.o p5namot-vdw.o p5namot-Bases.o p5namot-graphics.o p5namot-read_lib.o p5namot-spit.o p5namot-vrml.o p5namot-builder.o p5namot-scriptinter.o p5namot-read_parm.o p5namot-test.o p5namot-cmd_line.o p5namot-input.o p5namot-pairing.o p5namot-read_pdb.o p5namot-cmds_h.o p5namot-mathutil.o p5namot-input2.o p5namot-pcode.o p5namot-render.o p5namot-utility.o p5namot-pngwriter.o p5namot-strlcat.o p5namot-strlcpy.o p5namot-hub.o p5namot-mcg.o p5namot-x11stubs.o p5namot-guistubs.o p5namot-brokenlinks.o p5namot-snapshot.o p5namot-vectors.o p5namot-singlestrand.o p5namot-amberff.o p5namot-mc.o p5namot-read_amber7.o '  # Object files
     );
