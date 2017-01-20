# MB-RRL2 
### a library for model-based relational reinforcement learning


mbrrl2 is a fork of [libPRADA](http://userpage.fu-berlin.de/tlang/prada/)

If you have questions, please contact [Sebastian Höfer](mailto:science@sebastianhoefer.de).

#### Main features of this fork

- New object-oriented interface of the relational rule learner
- Task-sensitive learning as presented in paper [1]
- Incremental rule learning
- Several bugfixes

### Installation

Please type "make" in the root directory mbrrl2/.

Then, try out the demos in mbrrl2/test/.

In case of problems, please contact me.

to understand the main features of this library, please read the guide:

   doc/guide.pdf


### When using this library...

... please cite the following paper:

    [1] Sebastian Höfer, Tobias Lang, Oliver Brock
    Extracting Kinematic Background Knowledge from Interactions Using Task-Sensitive Relational Learning. 
    Proceedings of the IEEE International Conference on Robotics and Automation (ICRA), 
    2014, pp. 4342-4347.

Bibtex:

    @inproceedings{hoefer_extracting_2014,
      Title = {Extracting Kinematic Background Knowledge from Interactions Using Task-Sensitive Relational Learning},
      Author = {Sebastian Höfer and Tobias Lang and Oliver Brock},
      Booktitle = {Proceedings of the IEEE International Conference on Robotics and Automation (ICRA)},
      Pages = {4342-4347},
      Year = {2014},
      Doi = {10.1109/ICRA.2014.6907491},
      Location = {Hong Kong, China},
      Url = {http://www.robotics.tu-berlin.de/fileadmin/fg170/Publikationen_pdf/hoefer_extracting_2014.pdf},
    }


If you use the planner PRADA, please cite:

    [2] Planning with Noisy Probabilistic Relational Rules
    Journal of Artificial Intelligence Research
    2010, Volume 39, pages 1-49  

Bibtex:

    @article{lang-toussaint-10jair,
      author = {Tobias Lang and Marc Toussaint},
      title = {Planning with Noisy Probabilistic Relational Rules},
      journal = {Journal of Artificial Intelligence Research},
      year = {2010},
      volume = {39},
      pages = {1-49},
      pdfurl = {http://www.jair.org/media/3093/live-3093-5172-jair.pdf},
    }

If you use the relational rule learner, please cite the original publication presenting this learner:

    [3]  Pasula, Hanna M., Luke S. Zettlemoyer, and Leslie Pack Kaelbling. 
    Learning symbolic models of stochastic domains. 
    Journal of Artificial Intelligence Research 29 (2007): 309-352.

Bibtex:

    @article{pasula2007learning,
      title={Learning symbolic models of stochastic domains},
      author={Pasula, Hanna M and Zettlemoyer, Luke S and Kaelbling, Leslie Pack},
      journal={Journal of Artificial Intelligence Research},
      volume={29},
      pages={309--352},
      year={2007}
    }
    
### GPL licence statement:

mbrrl2 is a fork of libPRADA by Sebastian Höfer <science@sebastianhoefer.de>.

libPRADA Copyright 2008-2012 Tobias Lang
email: tobias.lang@fu-berlin.de

This file is part of mbrrl2/libPRADA.

mbrrl2/libPRADA is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

mbrrl2 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with mbrrl2.  If not, see <http://www.gnu.org/licenses/>

