from snakemake.shell import shell

db_user = "genome"
db_host = "genome-mysql.soe.ucsc.edu"
db_port = 3306
db_database = "hg38"

columns = snakemake.params.get("columns", "*")
table = snakemake.params.get("table")
where_regexp = snakemake.params.get("where_regexp", None)
order_by = snakemake.params.get("order_by", None)

extra = snakemake.params.get("extra", "")

sql_statement = "SELECT {} FROM {}".format(columns, table)

if where_regexp is not None:
    sql_statement = " ".join(
        [sql_statement, "WHERE {} REGEXP '{}'".format(where_regexp[0], where_regexp[1])]
    )

if order_by is not None:
    sql_statement = " ".join([sql_statement, "ORDER BY {}".format(order_by)])

shell(
    "mysql "
    "--user={db_user} "
    "--host={db_host} "
    "-P {db_port} "
    "-D {db_database} "
    "{extra} "
    "-A -B "
    '-e "{sql_statement};" '
    "> {snakemake.output}"
)
