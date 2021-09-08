function load_database(tm_info) {
var table = new Tabulator("#tm-results-table", {
    data:tm_info, //assign data to table
columns:[
{title:"Id", field:"TEMPLATE_MATCH_ID", sorter:"number"},
{title:"Image", field:"IMAGE_NAME", sorter:"string"},
{title:"Template", field:"REFERENCE_NAME", sorter:"string"},
{title:"Matches", field:"NUM_MATCHES", sorter:"number"},
{title:"Mean", field:"AVG_SCORE", sorter:"number",formatter:"money",formatterParams:{

precision:2,}},
{title:"Max", field:"MAX_SCORE", sorter:"number",formatter:"money",formatterParams:{

precision:2,
}},
], selectable:1, rowSelected:function(row){
window.location = "templatematching:load_job_from_database";

}         });
}