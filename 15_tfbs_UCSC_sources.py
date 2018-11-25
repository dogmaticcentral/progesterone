#!/usr/bin/python3

#
# This file is part of Progesternoe pipeline.
#
# Progesterone pipeline  is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Progesterone pipeline is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Progesterone pipeline.  If not, see <https://www.gnu.org/licenses/>.
#

from utils.mysqldb import *

# finding exps by exp num in UCSC:
#  SET @row_number = 0;  select * from (SELECT (@row_number:=@row_number + 1) \
#  AS expNum, source, factor, cellType, treatment, lab from wgEncodeRegTfbsClusteredInputsV3) as t where factor='ESR1';
# (the wgEncodeRegTfbsClusteredInputsV3 table has no id nums, so we have to count the rows)
#########################################
def main():

	assembly = "hg19" # afaik this is the only assembly with ENCODE data
	if len(sys.argv) < 2:
		print  ("usage:   %s <TF name> " % sys.argv[0])
		exit()

	tf_name = sys.argv[1]
	db     = connect_to_mysql("/home/ivana/.ucsc_mysql_conf")
	cursor = db.cursor()
	switch_to_db(cursor, assembly) # human build name

	# our table du jour is wgEncodeRegTfbsClusteredV3;
	table = 'wgEncodeRegTfbsClusteredInputsV3'

	qry = "SET @row_number = 0"
	search_db(cursor,qry)
	subquery  = "SELECT (@row_number:=@row_number + 1) "
	subquery += "AS expNum, source, factor, cellType, treatment, lab from %s " % table
	qry = "select * from (%s) as t  where factor='%s'" % (subquery, tf_name)
	# columns: tableName, source, factor, antibody, cellType, treatment, lab
	ret = search_db(cursor,qry)

	if ret==None:
		print ("No ret for %s" %qry)
		exit()
	if not ret or  isinstance([0][0],str) and 'Error'in ret[0][0]:
		search_db(cursor, qry, verbose=True)
		exit()
	print("\t".join(['tableName', 'source', 'factor', 'antibody', 'cellType', 'treatment', 'lab']))
	print("\n".join(["\t".join([str(field) for field in row]) for row in ret]))

	cursor.close()
	db.close()


	return True


#########################################
########################################
if __name__ == '__main__':
	main()

