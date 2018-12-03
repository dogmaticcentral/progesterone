#
# This file is part of Progesterone pipeline.
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

import MySQLdb
import sys

########
def connect_to_mysql (conf_file):
	try:
		mysql_conn_handle = MySQLdb.connect(read_default_file=conf_file)
	except  MySQLdb.Error as e:
		print("Error connecting to mysql (%s) " % (e.args[1]))
		sys.exit(1)
	return mysql_conn_handle

########
def switch_to_db(cursor, db_name):
	qry = "use %s" % db_name
	rows = search_db(cursor, qry, verbose=False)
	if (rows):
		print(rows)
		return False
	return True


#######
def search_db(cursor, qry, verbose=False):
	try:
		cursor.execute(qry)
	except MySQLdb.Error as e:
		if verbose:
			print("Error running cursor.execute() for  qry:\n%s\n%s" % (qry, e.args[1]))
		return [["Error"], e.args]

	try:
		rows = cursor.fetchall()
	except MySQLdb.Error as e:
		if verbose:
			print("Error running cursor.fetchall() for  qry:\n%s\n%s" % (qry, e.args[1]))
		return [["Error"], e.args]

	if len(rows) == 0:
		if verbose:
			print("No return for query:\n%s" % qry)
		return False

	return rows


########
def hard_check (db,cursor, ret, qry):
	if not ret or (type(ret[0][0])==str and 'Error' in ret[0][0]):
		search_db(cursor,qry, verbose=True)
		cursor.close()
		db.close()
		exit()


########
def permissive_check (db,cursor, ret, qry):
	if ret and (type(ret[0][0])==str and 'Error' in ret[0][0]):
		search_db(cursor,qry, verbose=True)
		cursor.close()
		db.close()
		exit()

########
def val2mysqlval(value):
	if  value is None:
		return  "null "
	elif type(value) is str:
		return "\'%s\'" % value
	return "{}".format(value)


########
def store_without_checking(cursor, table, fields, verbose=False):
	qry = "insert into %s " % table
	qry += "("
	qry += ",".join(fields.keys())
	qry += ")"

	qry += " values "
	qry += "("
	qry += ",".join([val2mysqlval(v) for v in fields.values()])
	qry += ")"

	rows  = search_db (cursor, qry, verbose)
	if verbose: print("qry:",qry,"\n", "rows:", rows)

	if rows:
		rows   = search_db (cursor, qry, verbose=True)
		print(rows)
		return -1

	rows = search_db (cursor, "select last_insert_id()" )
	try:
		row_id = int(rows[0][0])
	except:
		row_id = -1
	return row_id

########
def store_or_update (cursor, table, fixed_fields, update_fields, verbose=False, primary_key='id'):

	conditions = " and ".join(["{}={}".format(k,val2mysqlval(v)) for k,v in fixed_fields.items()])

	# check if the row exists
	qry = "select %s from %s  where %s "  % (primary_key, table, conditions)
	rows   = search_db (cursor, qry, verbose)
	exists = rows and (type(rows[0][0]) is int)

	row_id = -1
	if exists: row_id = rows[0][0]
	if verbose: print("\n".join(["", qry, "exists? {}".format(exists), str(row_id)]))
	if exists and not update_fields: return row_id

	if exists: # if it exists, update
		if verbose: print("exists; updating")
		qry  = "update %s set " % table
		qry += ",".join(["{}={}".format(k,val2mysqlval(v)) for k,v in update_fields.items()])
		qry += " where %s " % conditions

	else: # if not, make a new one
		if verbose: print("does not exist; making new one")
		qry  = "insert into %s " % table
		keys = list(fixed_fields.keys())
		vals = list(fixed_fields.values())
		if update_fields:
			keys += list(update_fields.keys())
			vals += list(update_fields.values())
		qry += "(" + ",".join(keys) + ")"
		qry += " values "
		qry += "(" + ",".join([val2mysqlval(v) for v in vals]) + ")"

	rows   = search_db (cursor, qry, verbose)

	if verbose: print("qry:",qry,"\n", "rows:", rows)
	# if there is a return, it is an error msg
	if rows:
		rows   = search_db (cursor, qry, verbose=True)
		print(rows[0])
		return -1

	if row_id==-1:
		rows = search_db (cursor, "select last_insert_id()" )
		try:
			row_id = int(rows[0][0])
		except:
			row_id = -1
	return row_id


#########################################
def store_xref (cursor, xtype, xid, bibtex=None, parent_id=None ):
	fixed_fields = {'xtype': xtype,'external_id': xid}
	if bibtex or parent_id:
		update_fields={}
		if bibtex: update_fields['bibtex'] = bibtex
		if parent_id: update_fields['parent_id'] = parent_id
	else:
		update_fields = None
	xref_id = store_or_update (cursor, 'xrefs', fixed_fields, update_fields)
	if xref_id<0:
		print ("Error storing", xtype, xid)
		exit()
	return xref_id


#########################################
def get_xref_id (db, cursor, external_exp_id):
	qry = "select id from xrefs where external_id='%s'" % external_exp_id
	ret = search_db(cursor,qry)
	hard_check (db, cursor, ret, qry)
	return int(ret[0][0])

#########################################
def get_gene_region_id (db, cursor, gene_name, assembly):
	qry = "select r.id  from regions as r, genes as g "
	qry += "where g.name='%s' and g.region_id=r.id and r.assembly='%s' " % (gene_name,assembly)
	ret = search_db(cursor,qry)
	hard_check (db,cursor, ret, qry)
	return ret[0][0]

#########################################
def get_region_coords (db, cursor, region_id):
	qry  = "select chromosome, rfrom, rto, strand from regions "
	qry += "where id=%d" % region_id
	ret = search_db(cursor,qry)
	hard_check (db, cursor, ret, qry)
	return ret[0]


#########################################
def get_gene_coords (db, cursor, gene_name, assembly):
	qry = "select r.chromosome, r.rfrom, r.rto, r.strand from regions as r, genes as g "
	qry += "where g.name='%s' and g.region_id=r.id and r.assembly='%s' " % (gene_name,assembly)
	ret = search_db(cursor,qry)
	hard_check (db,cursor, ret, qry)
	# there might be multiple returns, corresponding to different splices
	[chromosome, min_start, max_end, strand] = ret[0]
	for row in ret:
		[chromosome, start, end, strand] = row
		min_start = start if min_start>start else min_start
		max_end = end if max_end<end else max_end
	return [chromosome, strand, min_start, max_end]


#########################################
def get_tad_region(db, cursor, exp_file_xref_id, chromosome, min_start, max_end):

	# TODO: what if the gene is stradling the TAD region?
	qry  = "select rfrom, rto from regions where xref_id=%d " % exp_file_xref_id
	qry += "and chromosome='%s' and rfrom<%d and %d<rto" % (chromosome, min_start, min_start)
	ret = search_db(cursor,qry)
	hard_check (db,cursor, ret, qry)
	return ret[0]


########################################
def get_all_tads(db, cursor, exp_file_xref_id, chromosome):
	qry  = "select rfrom, rto from regions where xref_id=%d " % exp_file_xref_id
	qry += "and chromosome='%s' " % chromosome
	ret = search_db(cursor,qry)
	hard_check (db,cursor, ret, qry)
	return ret

########################################
def get_binding_regions(db, cursor, assembly, chromosome, tf_name, return_binding_site_id=False):

	qry   = "select "
	if return_binding_site_id: qry  += "b.id, "
	qry  += "r.rfrom, r.rto from regions as r, binding_sites as b "
	qry  += "where b.tf_name='%s' " % tf_name
	qry  += "and b.region_id = r.id "
	qry  += "and r.assembly='%s' and r.chromosome='%s' " % (assembly, chromosome)
	ret = search_db(cursor,qry)
	hard_check (db,cursor, ret, qry)
	return ret


########################################
def get_binding_regions_in_interval(db, cursor, assembly, chromosome, interval_start, interval_end, tf_name, return_binding_site_id=False):
	qry   = "select "
	if return_binding_site_id: qry  += "b.id, "
	qry  += "r.rfrom, r.rto from regions as r, binding_sites as b "
	qry  += "where b.tf_name='%s' " % tf_name
	qry  += "and b.region_id = r.id "
	qry  += "and r.assembly='%s' and r.chromosome='%s' " % (assembly, chromosome)
	qry  += "and r.rfrom>=%d and r.rto<=%d " % (interval_start, interval_end)
	ret = search_db(cursor,qry)
	hard_check (db,cursor, ret, qry)
	return ret


########################################
def get_motifs_in_binding_site(db, cursor, binding_site_id):
	qry = "select motif_id from binding_site2motif "
	qry += "where binding_site_id='%d' " % binding_site_id
	ret = search_db(cursor,qry)
	permissive_check(db,cursor, ret, qry)
	return [r[0] for r in ret] if ret else []

########################################
def assembly2species_common(cursor,assembly):
	ret = search_db(cursor,"select common_name from assemblies where  assembly='%s'"% assembly)
	if not ret or type(ret[0][0])!=str or 'Error' in ret[0][0]:
		return ""
	return ret[0][0].strip().lower()
