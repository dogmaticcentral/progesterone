
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
			print("Error running cursor.execute() for  qry: %s: %s " % (qry, e.args[1]))
		return ["ERROR: " + e.args[1]]

	try:
		rows = cursor.fetchall()
	except MySQLdb.Error as e:
		if verbose:
			print("Error running cursor.fetchall() for  qry: %s: %s " % (qry, e.args[1]))
		return ["ERROR: " + e.args[1]]

	if len(rows) == 0:
		if verbose:
			print("No return for query %s" % qry)
		return False

	return rows


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
		return False

	return True

########
def store_or_update (cursor, table, fixed_fields, update_fields, verbose=False, primary_key='id'):

	conditions = "and ".join(["{}={}".format(k,val2mysqlval(v)) for k,v in fixed_fields.items()] )

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
		qry += ",".join(["{}={}".format(k,val2mysqlval(v)) for k,v in update_fields.items()] )
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
def store_xref(cursor, xtype, xid, bibtex=None ):
	fields = {'xtype': xtype,'external_id': xid}
	update_fields={'bibtex':bibtex} if bibtex else None
	xref_id = store_or_update (cursor, 'xrefs', fields, update_fields)
	if xref_id<0:
		print ("Error storing", xtype, xid)
		exit()
	return xref_id
