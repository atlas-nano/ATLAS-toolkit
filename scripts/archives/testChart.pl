#!/usr/bin/perl -w
BEGIN {
    unshift @INC, "/ul/tpascal/scripts/Packages";
}

use Chart::Graph::Gnuplot qw(gnuplot);

     gnuplot({"title" => "Examples of Errorbars",
              "xrange" => "[:11]",
              "yrange" => "[:45]",
              "output file" => "gnuplot2.gif",
              "output type" => "gif",
             },
             # dataset 1
             [{"title" => "yerrorbars",
               "style" => "yerrorbars",
               "using" => "1:2:3:4",
               "type" => "columns"},
              [ 1, 2, 3, 4, 5, 6 ], # x
              [ 5, 7, 12, 19, 28, 39 ], # y
              [ 3, 5, 10, 17, 26, 38 ], # ylow
              [ 6, 8, 13, 20, 30, 40 ] ], # yhigh
             # dataset 2
             [{"title" => "xerrorbars",
               "style" => "xerrorbars",
               "using" => "1:2:3:4",
               "type" => "columns"},
              [ 4, 5, 6, 7, 8, 9 ], # x
              [ 1, 4, 5, 6, 7, 10 ], # y
              [ 3.3, 4.4, 5.5, 6.6, 7.7, 8.8 ], # xlow
              [ 4.1, 5.2, 6.1, 7.3, 8.1, 10 ] ], # xhigh
             # dataset 3
             [{"title" => "xyerrorbars",
               "style" => "xyerrorbars",
               "using" => "1:2:3:4:5:6",
               "type" => "columns"},
              [ 1.5, 2.5, 3.5, 4.5, 5.5, 6.5 ], # x
              [ 2, 3.5, 7.0, 14, 15, 20 ], # y
              [ 0.9, 1.9, 2.8, 3.7, 4.9, 5.8 ], # xlow
              [ 1.6, 2.7, 3.7, 4.8, 5.6, 6.7 ], # xhigh
              [ 1, 2, 3, 5, 7, 8 ], # ylow
              [ 5, 7, 10, 17, 18, 24 ] ], # yhigh
             # dataset 4
             [{"title" => "xerrorbars w/ xdelta",
               "style" => "xerrorbars",
               "using" => "1:2:3",
               "type" => "columns"},
              [ 4, 5, 6, 7, 8, 9 ], # x
              [ 2.5, 5.5, 6.5, 7.5, 8.6, 11.7 ], # y
              [ .2, .2, .1, .1, .3, .3 ] ], # xdelta
             # dataset 5
             [{"title" => "yerrorbars w/ ydelta",
               "style" => "yerrorbars",
               "using" => "1:2:3",
               "type" => "columns"},
              [ .7, 1.7, 2.7, 3.7, 4.7, 5.7 ], # x
              [ 10, 15, 20, 25, 30, 35 ], # y
              [ .8, 1.2, 1.1, 2.1, 1.3, 3.3 ] ], # ydelta
             # dataset 6
             [{"title" => "dummy data",
               "type" => "matrix"},
              [ [1,1] ]],
             # dataset 7
             [{"title" => "xyerrorbars w/ xydelta",
               "style" => "xyerrorbars",
               "using" => "1:2:3:4",
               "type" => "columns"},
               [ 7.5, 8.0, 8.5, 9.0, 9.5, 10.0 ], # x
               [ 30, 27, 25, 23, 27, 33 ], # y
               [ .2, .1, .3, .6, .4, .3 ], # xdelta
              [ .8, .7, .3, .6, 1.0, .3 ] ], # ydelta
           );
gnuplot({'title' => 'Corporate stock values for a major computer maker',
           'x-axis label' => 'Month and Year',
           'y-axis label' => 'Stock price',
           'output type' => 'png',
           'output file' => 'gnuplot3.png',
           # Setting date/time specific options.
           'xdata' => 'time',
           'timefmt' => '%m/%d/%Y',
           'format' => ['x', '%m/%d/%Y'],
           # Set output range - note quoting of date string
           'xrange' => '["06/01/2000":"08/01/2001"]',
           'extra_opts' => join("\n", 'set grid', 'set timestamp'),
          },
          # Data for when stock opened
          [{'title' => 'open',
            'type' => 'matrix',
            'style' => 'lines',
           },
           [
            ['06/01/2000',  '81.75'],
            ['07/01/2000', '52.125'],
            ['08/01/2000', '50.3125'],
            ['09/01/2000', '61.3125'],
            ['10/01/2000', '26.6875'],
            ['11/01/2000', '19.4375'],
            ['12/01/2000', '17'],
            ['01/01/2001', '14.875'],
            ['02/01/2001', '20.6875'],
            ['03/01/2001', '17.8125'],
            ['04/01/2001', '22.09'],
            ['05/01/2001', '25.41'],
            ['06/01/2001', '20.13'],
            ['07/01/2001', '23.64'],
            ['08/01/2001', '19.01'],
           ]
          ],

          # Data for stock high
          [{'title' => 'high',
            'type' => 'matrix',
            'style' => 'lines',
           },
           [
            ['06/01/2000', '103.9375'],
            ['07/01/2000', '60.625'],
            ['08/01/2000', '61.50'],
            ['09/01/2000', '64.125'],
            ['10/01/2000', '26.75'],
            ['11/01/2000', '23'],
            ['12/01/2000', '17.50'],
            ['01/01/2001', '22.50'],
            ['02/01/2001', '21.9375'],
            ['03/01/2001', '23.75'],
            ['04/01/2001', '27.12'],
            ['05/01/2001', '26.70'],
            ['06/01/2001', '25.10'],
            ['07/01/2001', '25.22'],
            ['08/01/2001', '19.90'],
           ]
          ],

          # Data for stock close
          [{'title' => 'close',
            'type' => 'matrix',
            'style' => 'lines',
           },
           [

            ['06/01/2000', '52.375'],
            ['07/01/2000', '50.8125'],
            ['08/01/2000', '60.9375'],
            ['09/01/2000', '25.75'],
            ['10/01/2000', '19.5625'],
            ['11/01/2000', '16.50'],
            ['12/01/2000', '14.875'],
            ['01/01/2001', '21.625'],
            ['02/01/2001', '18.25'],
            ['03/01/2001', '22.07'],
            ['04/01/2001', '25.49'],
            ['05/01/2001', '19.95'],
            ['06/01/2001', '23.25'],
            ['07/01/2001', '18.79'],
            ['08/01/2001', '18.55'],
           ]
          ]
);
