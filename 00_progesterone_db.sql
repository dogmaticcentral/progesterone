--
--
-- This file is part of Progesterone pipeline.
--
-- Progesterone pipeline  is free software: you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation, either version 3 of the License, or
-- (at your option) any later version.
--
-- Progesterone pipeline is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
--
-- You should have received a copy of the GNU General Public License
-- along with Progesterone pipeline.  If not, see <https://www.gnu.org/licenses/>.
--
--
--

CREATE DATABASE progesterone;
USE progesterone;

CREATE TABLE xrefs (
  id int  NOT NULL AUTO_INCREMENT PRIMARY KEY,
  xtype ENUM ('pubmed','geo','encode', 'ucsc') NOT NULL,
  external_id varchar(255) NOT NULL,
  parent_id int DEFAULT NULL,
  bibtex text DEFAULT NULL
) ENGINE=InnoDB;

CREATE TABLE regions (
  id int  NOT NULL AUTO_INCREMENT PRIMARY KEY,
  species  varchar(50) NOT NULL,
  chromosome varchar(50) NOT NULL,
  assembly varchar(50) NOT NULL,
  rfrom bigint  NOT NULL,
  rto bigint  NOT NULL,
  strand varchar(1) DEFAULT NULL,
  rtype ENUM ('chromosome','gene', 'tad', 'interacting', 'motif', 'chipseq', 'atacseq') NOT NULL,
  xref_id int DEFAULT NULL,
  FOREIGN KEY fk_xref(xref_id) REFERENCES xrefs(id) ON UPDATE CASCADE ON DELETE SET NULL
) ENGINE=InnoDB;

-- have multiple entries if there are multiple regions referring to the gene - species, assemblies
CREATE TABLE  genes (
  id int  NOT NULL AUTO_INCREMENT PRIMARY KEY,
  name varchar(50) NOT NULL,
  region_id  int DEFAULT NULL,
  FOREIGN KEY fk_region(region_id) REFERENCES regions(id) ON UPDATE CASCADE ON DELETE CASCADE
) ENGINE=InnoDB;

CREATE TABLE interactions (
  id int  NOT NULL AUTO_INCREMENT PRIMARY KEY,
  region_id1 int NOT NULL,
  region_id2 int NOT NULL,
  ivalue int NOT NULL,
  xref_id int DEFAULT NULL,
  FOREIGN KEY fk_region1(region_id1) REFERENCES regions(id) ON UPDATE CASCADE ON DELETE CASCADE,
  FOREIGN KEY fk_region2(region_id2) REFERENCES regions(id) ON UPDATE CASCADE ON DELETE CASCADE,
  FOREIGN KEY fk_xref(xref_id) REFERENCES xrefs(id) ON UPDATE CASCADE ON DELETE CASCADE
) ENGINE=InnoDB;

-- motif ids and alignment use ";" as separator - their entries should match respectively
-- i.e the sequence in the alignment should be the same as in motifs, up to indels
CREATE TABLE alignments (
  id int  NOT NULL AUTO_INCREMENT PRIMARY KEY,
  motif_ids text NOT NULL,
  alignment text NOT NULL,
  xref_id int DEFAULT NULL,
  FOREIGN KEY fk_xref(xref_id) REFERENCES xrefs(id) ON UPDATE CASCADE ON DELETE SET NULL
) ENGINE=InnoDB;


CREATE TABLE motifs (
  id int  NOT NULL AUTO_INCREMENT PRIMARY KEY,
  region_id int NOT NULL,
  tf_name   varchar(50) NOT NULL,
  sequence  varchar(255) NOT NULL,
  consensus varchar(255) NOT NULL,
  score float DEFAULT NULL,
  xref_id      int NOT NULL,
  alignment_id int DEFAULT NULL,
  FOREIGN KEY fk_region(region_id) REFERENCES regions(id) ON UPDATE CASCADE ON DELETE CASCADE,
  FOREIGN KEY fk_xref(xref_id) REFERENCES xrefs(id) ON UPDATE CASCADE ON DELETE CASCADE,
  FOREIGN KEY fk_almt(alignment_id) REFERENCES alignments(id) ON UPDATE CASCADE ON DELETE SET NULL
) ENGINE=InnoDB;


CREATE TABLE  binding_sites (
  id int  NOT NULL AUTO_INCREMENT PRIMARY KEY,
  tf_name   varchar(50) NOT NULL,
  motif_ids text DEFAULT NULL,
  region_id int NOT NULL,
  xref_id int DEFAULT NULL,
  FOREIGN KEY fk_region(chipseq_region_id) REFERENCES regions(id) ON UPDATE CASCADE ON DELETE CASCADE,
  FOREIGN KEY fk_xref(xref_id) REFERENCES xrefs(id) ON UPDATE CASCADE ON DELETE CASCADE
) ENGINE=InnoDB;

CREATE TABLE  binding_site2motif (
  id int  NOT NULL AUTO_INCREMENT PRIMARY KEY,
  binding_site_id int  NOT NULL,
  motif_id int  NOT NULL,
  FOREIGN KEY fk_bs(binding_site_id) REFERENCES binding_sites(id) ON UPDATE CASCADE ON DELETE CASCADE,
  FOREIGN KEY fk_mf(motif_id) REFERENCES motifs(id) ON UPDATE CASCADE ON DELETE CASCADE
) ENGINE=InnoDB;


