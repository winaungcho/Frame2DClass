

<html>
    <head>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>Frame2D Analysis</title>
<style>
.fixedheadertable {
  width: 100%;
  max-height: 60%;
  overflow: scroll;
}

table {
  position: relative;
  border: 1px solid #ddd;
  border-collapse: collapse;
    text-decoration:none;
}

td, th {
  white-space: nowrap;
  border: 1px solid #ddd;
  padding: 5px;
  
}

th {
  background-color: #eee;
  position: -webkit-sticky;
  position: sticky;
  top: -1px;
  z-index: 2;
  text-align: center;
}
th:first-of-type {
  left: 0;
  z-index: 3;
}

tbody tr td:first-of-type {
  background-color: #eee;
  position: -webkit-sticky;
  position: sticky;
  left: -1px;
  text-align: left;
}
tr:nth-child(even) {
  background-color: #fafafa;
}
.hover {
  background: yellow;
}
</style>
    </head>
    <body>
<?php
include("frame2dclass.php");

$pframe = new Frame2D();
$pframe->Init();
//echo " init";
echo $pframe->Process();
$pframe->saveCSV("frame2d.csv");
$pframe->readGenCSV("frame2d.csv");
echo $pframe->Process();
?>
  
    </body>
</html>
