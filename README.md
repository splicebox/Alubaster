## Alubaster
Detection of novel *Alu* exonization events from RNA-seq data

Described in:
  Florea L, Payer L, Antonescu C, Yang G and Burns K. (2021) Detection of Alu exonization events in human frontal cortex from RNA-seq data, *Submitted*.

Supplementary data described in the article can be found in the '[alubaster_paper](https://github.com/splicebox/ALubaster/tree/master/alubaster_paper)' directory.

```
Copyright (C) 2018-2021 and GNU GPL v3.0 by Liliana Florea, Lindsay Payer
```

### Synopsis
Alubaster identifies candidate gene loci by locating mapped reads (‘anchors’) whose *Alu*-containing mates could not be found in the genome. It then applies two types of filters to select a more accurate subset of loci. The first, a *signal filter*, identifies read evidence for an Alu exonization event, by searching the mate’s sequence against a concatenation of the neighboring exons’ and *Alu* sequences while concomitantly ruling out false positive matches that could have resulted from local or more distant *Alu* elements in the genome. The second, a *context filter*, evaluates the likelihood of an event based on the strength of the signal versus the local context, in particular repeat content and proportion of signal-to-‘context’ matches. 

### Version
The software is alpha version. A user-friendly version will be available soon.

### Terms of use
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received (LICENSE.txt) a copy of the GNU General Public License v3.0 along with this program; if not, you can obtain one from http://www.gnu.org/licenses/gpl.txt or by writing to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
