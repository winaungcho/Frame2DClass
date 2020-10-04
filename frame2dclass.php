<?php
/******
 * Frame2DClass
 *
 * Class Inheritance
 * [FEMSolver]
 *  -> [Frame2D]
 * [FEMSolver] is base class and includes several matrix operation for the standard FEM solutions.
 * [Frame2D] is a class for the FEM solution process and include data structure of 2 dimensional rigid jointed Frame.
 * [Frame2D] use [CSV] class for handling csv file to store and retrieve associated arrays;
 * Solution process run for the loaded truss to analyse deformations, reactions and element forces.
 * Multiple load cases will be solved simultaneously.
 * Html result tables are generated during the process.
 * Model of 2d Frame can be generated within class by assigning values to variables.
 * Or model can be created by loading CSV file.
 * Model csv file is very simple comma separated text file in-wich the properties of FEM element, boundary conditions and loads are written.
 * 
 * This class is free for the educational use as long as maintain this header together with this class.
 * Author: Win Aung Cho
 * Contact winaungcho@gmail.com
 * version 1.0
 * Date: 30-9-2020
 *
 ******/

class FEMSolver extends CSV
{
    function lambada($dx, $dy)
    {
        //Rotation Matrix 2d, x, y, θz
        $ELL = sqrt(pow(($dy) , 2) + pow(($dx) , 2));
        $dcos[0][0] = $dx / $ELL;
        $dcos[0][1] = $dy / $ELL;
        $dcos[0][2] = 0;
        $dcos[1][0] = - $dcos[0][1];
        $dcos[1][1] = $dcos[0][0];
        $dcos[1][2] = 0;
        $dcos[2][0] = 0;
        $dcos[2][1] = 0;
        $dcos[2][2] = 1.0;
        for ($k = 0;$k < 6;$k++)
        {
            for ($j = 0;$j < 6;$j++)
            {
                $alambda[$k][$j] = 0;
            }
        }
        for ($k = 0;$k < 2;$k++)
        {
            $ik = 3 * $k;
            for ($i = 0;$i < 3;$i++)
            {
                for ($j = 0;$j < 3;$j++)
                {
                    $alambda[$i + $ik][$j + $ik] = $dcos[$i][$j];
                }
            }
        }
        return $alambda;
    }
    function BmCol1(&$Km, $EA, $EI, $dx, $dy)
    {
        // Klocal Beam2d FxI, FyI, MzI, FxJ, FyJ, MzJ
        $ELL = sqrt(pow(($dy) , 2) + pow(($dx) , 2));
        for ($i = 0;$i < 6;$i++)
        {
            for ($j = $i;$j < 6;$j++)
            {
                $Km[$i][$j] = 0;
            }
        }
        $Km[0][0] = $EA / $ELL;
        $Km[0][3] = - $EA / $ELL;
        $Km[3][3] = $EA / $ELL;
        $Km[1][1] = 12.0 * $EI / $ELL / $ELL / $ELL;
        $Km[1][2] = 6.0 * $EI / $ELL / $ELL;
        $Km[1][4] = - $Km[1][1];
        $Km[1][5] = $Km[1][2];
        $Km[2][2] = 4.0 * $EI / $ELL;
        $Km[2][4] = - 6.0 * $EI / $ELL / $ELL;
        $Km[2][5] = 2.0 * $EI / $ELL;
        $Km[4][4] = 12.0 * $EI / $ELL / $ELL / $ELL;
        $Km[4][5] = - 6.0 * $EI / $ELL / $ELL;
        $Km[5][5] = 4.0 * $EI / $ELL;
        for ($i = 0;$i < 6;$i++)
        {
            for ($j = $i;$j < 6;$j++)
            {
                if (isset($Km[$i][$j])) $Km[$j][$i] = $Km[$i][$j];
            }
        }
    }
    function BmCol2(&$Km, $EA, $EI, $dx, $dy)
    {
        //Kglobal Beam2d FxI, FyI, MzI, FxJ, FyJ, MzJ
        //float $ELL, $C, $S, $E1, $E2, $E3, $E4;
        $ELL = sqrt(pow(($dy) , 2) + pow(($dx) , 2));
        $C = ($dx) / $ELL;
        $S = ($dy) / $ELL;
        $E1 = $EA / $ELL;
        $E2 = 12.0 * $EI / $ELL / $ELL / $ELL;
        $E3 = $EI / $ELL;
        $E4 = 6.0 * $EI / $ELL / $ELL;
        $Km[0][0] = $C * $C * $E1 + $S * $S * $E2;
        $Km[3][3] = $Km[0][0];
        $Km[0][1] = $S * $C * ($E1 - $E2);
        $Km[1][0] = $Km[0][1];
        $Km[3][4] = $Km[0][1];
        $Km[4][3] = $Km[3][4];
        $Km[0][2] = - $S * $E4;
        $Km[2][0] = $Km[0][2];
        $Km[0][5] = $Km[0][2];
        $Km[5][0] = $Km[0][5];
        $Km[2][3] = $S * $E4;
        $Km[3][2] = $Km[2][3];
        $Km[3][5] = $Km[2][3];
        $Km[5][3] = $Km[3][5];
        $Km[0][3] = - $Km[0][0];
        $Km[3][0] = $Km[0][3];
        $Km[0][4] = $S * $C * ($E2 - $E1);
        $Km[4][0] = $Km[0][4];
        $Km[1][3] = $Km[0][4];
        $Km[3][1] = $Km[1][3];
        $Km[1][1] = $S * $S * $E1 + $C * $C * $E2;
        $Km[4][4] = $Km[1][1];
        $Km[1][4] = - $Km[1][1];
        $Km[4][1] = $Km[1][4];
        $Km[1][2] = $C * $E4;
        $Km[2][1] = $Km[1][2];
        $Km[1][5] = $Km[1][2];
        $Km[5][1] = $Km[1][5];
        $Km[2][2] = 4.0 * $E3;
        $Km[5][5] = $Km[2][2];
        $Km[2][4] = - $C * $E4;
        $Km[4][2] = $Km[2][4];
        $Km[4][5] = $Km[2][4];
        $Km[5][4] = $Km[4][5];
        $Km[2][5] = 2.0 * $E3;
        $Km[5][2] = $Km[2][5];
        return 1;
    }
    function assemble($node, $nbw, $Km, &$s)
    {
        //Assemble Stiffness into Global S
        $ndn = 3;
        $nen = 2;
        for ($ii = 0;$ii < $nen;$ii++)
        {
            $nrt = $ndn * ($node[$ii]);
            for ($it = 0;$it < $ndn;$it++)
            {
                $nr = $nrt + $it;
                $i = $ndn * $ii + $it;
                for ($jj = 0;$jj < $nen;$jj++)
                {
                    $nct = $ndn * ($node[$jj]);
                    for ($jt = 0;$jt < $ndn;$jt++)
                    {
                        $j = $ndn * $jj + $jt;
                        $nc = $nct + $jt - $nr;
                        if ($nc >= 0)
                        {
                            if (!isset($s[$nbw * $nr + $nc])) $s[$nbw * $nr + $nc] = 0;
                            $s[$nbw * $nr + $nc] += $Km[$i][$j];
                        }
                    }
                }
            }
        }
    }
    function modifydispbond($pred, &$s, $nbw, $cnst, &$f)
    {
        //modify for pre-defined displacement
        list($nd, $u, $nu) = $this->getdispbond($pred);
        // u = displacement, nu = index of dof, nd = no. of boundary cond.
        for ($i = 0;$i < $nd;$i++)
        {
            $k = $nu[$i];
            $s[($k) * $nbw] = $s[($k) * $nbw] + $cnst;
            if (!isset($f[$k])) $f[$k] = 0;
            $f[$k] = $f[$k] + $cnst * $u[$i];
        }
    }
    function getdispbond($pred)
    {
        // u = displacement, nu = index of dof, nd = no. of boundary cond.
        $nd = 0;
        $ndn = 3;
        $n = count($pred);
        for ($i = 0;$i < $n;$i++)
        {
            $nn = $pred[$i]["n"];
            if ($pred[$i]["ux"] == 0)
            {
                $u[$nd] = 0;
                $nu[$nd] = $ndn * $nn;
                $nd++;
            }
            if ($pred[$i]["uy"] == 0)
            {
                $u[$nd] = 0;
                $nu[$nd] = $ndn * $nn + 1;
                $nd++;
            }
            if ($pred[$i]["qz"] == 0)
            {
                $u[$nd] = 0;
                $nu[$nd] = $ndn * $nn + 2;
                $nd++;
            }
        }
        return array(
            $nd,
            $u,
            $nu
        );
    }
    function Reaction($pred, $cnst, $f)
    {
        //retrieve reaction
        $nd = 0;
        $ndn = 3;
        $n = count($pred);
        $reaction = array();
        for ($i = 0;$i < $n;$i++)
        {
            $nn = $pred[$i]["n"];
            $reaction[$i]["nn"] = $nn;
            if ($pred[$i]["ux"] == 0)
            {
                $k = $ndn * $nn;
                $reaction[$i]["ux"] = $cnst * (0 - $f[$k]);
            }
            if ($pred[$i]["uy"] == 0)
            {
                $k = $ndn * $nn + 1;
                $reaction[$i]["uy"] = $cnst * (0 - $f[$k]);
            }
            if ($pred[$i]["qz"] == 0)
            {
                $k = $ndn * $nn + 2;
                $reaction[$i]["qz"] = $cnst * (0 - $f[$k]);
            }
        }
        return $reaction;
    }
    function elementload($w, $dx, $dy)
    {
        // consistant distributed load in global
        $ed = array(
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0
        );
        $alambada = $this->lambada($dx, $dy);
        $ELL = sqrt(pow(($dy) , 2) + pow(($dx) , 2));
        $ed[0] = 0;
        $ed[3] = 0;
        $ed[1] = $w * $ELL / 2.0;
        $ed[4] = $ed[1];
        $ed[2] = $w * $ELL * $ELL / 12.0;
        $ed[5] = - $ed[2];
        $edp = array();
        for ($i = 0;$i < 6;$i++)
        {
            $edp[$i] = 0;
            for ($k = 0;$k < 6;$k++)
            {
                $edp[$i] = $edp[$i] + $alambada[$k][$i] * $ed[$k];
            }
        }
        return $edp;
    }
    function formelementload($node, $line, &$f, $loadcase, $ln)
    {
        //assemble consistant element load loadcase $ln
        $n = count($loadcase[$ln]["dists"]);
        for ($i = 0;$i < $n;$i++)
        {
            $ne = $loadcase[$ln]["dists"][$i]["n"];
            $node1 = $node[$line[$ne]["I"]];
            $node2 = $node[$line[$ne]["J"]];
            $dx = $node2["x"] - $node1["x"];
            $dy = $node2["y"] - $node1["y"];
            $w = $loadcase[$ln]["dists"][$i]["Wy"];
            $edp = $this->elementload($w, $dx, $dy);
            $i1 = 3 * $line[$ne]["I"];
            $i2 = 3 * $line[$ne]["J"];
            //echo $i1.":".$i2.":".$dx.":".$dy.":".$w;
            //print_r($edp);
            for ($j = 0;$j < 3;$j++)
            {
                if (!isset($f[$i1 + $j])) $f[$i1 + $j] = 0;
                if (!isset($f[$i2 + $j])) $f[$i2 + $j] = 0;
                $f[$i1 + $j] = $f[$i1 + $j] + $edp[$j];
                $f[$i2 + $j] = $f[$i2 + $j] + $edp[$j + 3];
            }
        };
    }
    function formf(&$f, $loadcase, $ln)
    {
        //assemble nodal load for loadcase $ln
        $ndn = 3;
        $nodekey = array("Fx","Fy","Mz");
        $n = count($loadcase[$ln]["nodes"]);
        for ($i = 0;$i < $n;$i++)
        {
            $nn = $loadcase[$ln]["nodes"][$i]["n"];
            for ($j = 0;$j < $ndn;$j++)
            {
                $nnq = $nn * $ndn + $j;
                if (!isset($f[$nnq])) $f[$nnq] = 0;
                $f[$nnq] = $f[$nnq] + $loadcase[$ln]["nodes"][$i][$nodekey[$j]];
            }
        };
    }
    function bansol($s, &$f, $nq, $nbw)
    {
        // linear solution of banded matrix
        //int $n1,$k,$nk,$i,$i1,$j,$j1,$kk;
        //float $c1;
        /* ----- band solver ----- */
        $n1 = $nq - 1;
        /* --- forward elimination --- */
        for ($k = 1;$k <= $n1;$k++)
        {
            $nk = $nq - $k + 1;
            if ($nk > $nbw) $nk = $nbw;
            for ($i = 2;$i <= $nk;$i++)
            {
                if (!isset($s[$nbw * ($k - 1) + $i - 1])) $s[$nbw * ($k - 1) + $i - 1] = 0;
                $c1 = $s[$nbw * ($k - 1) + $i - 1] / $s[$nbw * ($k - 1) ];
                $i1 = $k + $i - 1;
                for ($j = $i;$j <= $nk;$j++)
                {
                    $j1 = $j - $i + 1;
                    if (!isset($s[$nbw * ($k - 1) + $j - 1])) $s[$nbw * ($k - 1) + $j - 1] = 0;
                    if (!isset($s[$nbw * ($i1 - 1) + $j1 - 1])) $s[$nbw * ($i1 - 1) + $j1 - 1] = 0;
                    $s[$nbw * ($i1 - 1) + $j1 - 1] = $s[$nbw * ($i1 - 1) + $j1 - 1] - $c1 * $s[$nbw * ($k - 1) + $j - 1];
                }
                if (!isset($f[$i1 - 1])) $f[$i1 - 1] = 0;
                $f[$i1 - 1] = $f[$i1 - 1] - $c1 * $f[$k - 1];
            }
        }
        /* --- back-substitution --- */
        $f[$nq - 1] = $f[$nq - 1] / $s[$nbw * ($nq - 1) ];
        for ($kk = 1;$kk <= $n1;$kk++)
        {
            $k = $nq - $kk;
            $c1 = 1.0 / $s[$nbw * ($k - 1) ];
            $f[$k - 1] = $c1 * $f[$k - 1];
            $nk = $nq - $k + 1;
            if ($nk > $nbw) $nk = $nbw;
            for ($j = 2;$j <= $nk;$j++)
            {
                $f[$k - 1] = $f[$k - 1] - $c1 * $s[$nbw * ($k - 1) + $j - 1] * $f[$k + $j - 2];
            }
        }
    }
}
//echo " femsolve";
class Geometry
{
    function Dist($p1, $p2)
    {
        $dx = $p2["x"] - $p1["x"];
        $dy = $p2["y"] - $p1["y"];
        $L = sqrt($dx * $dx + $dy * $dy);
        return $L;
    }
    function Angle($p1, $p2)
    {
        $ang = 0;
        $dx = $p2["x"] - $p1["x"];
        $dy = $p2["y"] - $p1["y"];
        if ($dx == 0.0) if ($dy < 0.0) return -PI / 2.0;
        else if ($dy > 0.0) return PI / 2.0;
        $ang = atan(abs($dy / $dx));
        if (($dx < 0.0) && ($dy >= 0.0)) $ang = PI - $ang;
        else if (($dx < 0.0) && ($dy < 0.0)) $ang = PI + $ang;
        else if (($dx > 0.0) && ($dy < 0.0)) $ang = - $ang;
        return $ang;
    }
    function Polar($p, $dist, $ang)
    {
        $p["x"] += cos($ang) * $dist;
        $p["y"] += sin($ang) * $dist;
        return $p;
    }
    function DeflBeam($s, $l, $w1, $t1, $w2, $t2)
    {
        //var $N, $w;
        $N[0] = 1.0 - 3.0 * $s * $s + 2.0 * $s * $s * $s;
        $N[1] = $l * $s * (1.0 - 2.0 * $s + $s * $s);
        $N[2] = $s * $s * (3.0 - 2.0 * $s);
        $N[3] = $l * $s * $s * ($s - 1.0);
        $w = $N[0] * $w1 + $N[1] * $t1 + $N[2] * $w2 + $N[3] * $t2;
        return $w;
    }
    function GenDefl(&$p, $p1, $p2, $w1, $t1, $w2, $t2, $scale)
    {
        //var $w, $dL, $ang, $dx, $dy,$s;
        //var $i, $n=20;
        //var tmp;
        $n = 20;
        $dL = Dist($p1, $p2);
        $ang = Angle($p1, $p2);
        $dx = $p2["x"] - $p1["x"];
        $dy = $p2["y"] - $p1["y"];
        for ($i = 0;$i <= $n;$i++) $p[$i] = array(
            "x" => $p1["x"] + $dx / $n * ($i) ,
            "y" => $p1["y"] + $dy / $n * ($i)
        );
        if ($scale == 0) $scale = 200;
        for ($i = 0;$i <= $n;$i++)
        {
            $s = $i / 20.0;
            $w = DeflBeam($s, $dL, $w1, $t1, $w2, $t2);
            $p[$i] = Polar($p[$i], $w * $scale, $ang + PI / 2.0);
        }
        return $n + 1;
    }
}
class CSV {
    function __constructor(){
        
    }
    function saveArr2CSV($f, $data, $prefix){
        $style=1;
        // Add the keys as the column headers
        
        if ($style===1){
            fwrite($f, $prefix.",".count($data)."\r\n");
            fputcsv($f, array_keys($data[0]));
        }
        // Loop over the array and passing in the values only.
        foreach ($data as $row)
        {
            if ($style !== 1)
                fwrite($f, "$prefix,");
            fputcsv($f, $row);
        }
    }
    function readCSV($fname)
    {
        // open for reading
        $handle = fopen($fname, "r");
        if ($handle)
        {
            while (($head = fgetcsv($handle, 1000, ",")) !== false)
            {
                $nh = count($head);
                $n = (int) $head[$nh-1];
                //print_r($head)."<br/>";
                $name = $head[0];
                $arr = &$this->$name;
                if ($nh > 2)
                    $arr = &$arr[$head[1]];
                if ($nh > 3)
                    $arr = &$arr[$head[2]];
                if ($nh > 4)
                    $arr = &$arr[$head[3]];
                
                $data = fgetcsv($handle, 1000, ",");
                { // extract header data
                    $keys = $data; // save as keys
                    //print_r($keys);
                }
                for ($i = 0;$i < $n;$i++)
                {
                    $data = fgetcsv($handle, 1000, ",");
                    $arr[$i] = array_combine($keys, $data);
                    //print_r($arr[$i]);
                }
            }
            fclose($handle);
            //print_r($this->node);
        }
        else
        {
            // error opening the file.
            
        }
        
    }
}
class Frame2D extends FEMSolver
{
    var $node = array();
    var $boundary = array();
    var $line = array();
    var $Mat = array();
    var $Sec = array();
    var $nodalmass = array();
    var $loadcase = array();
    var $processfinished = false;
    var $linecolor = "#0000ff";
    var $textcolor = "#00ff00";
    var $resulthtml = "";
    // result
    var $react;
    function __constructor(){
        
    }
    function Init()
    {
        $this->node[] = array(
            "x" => 0.0,
            "y" => 0.0,
            "dof" => "000",
            "dx" => 0.00,
            "dy" => 0.00,
            "qz" => 0.0
        );
        $this->node[] = array(
            "x" => 0.0,
            "y" => 120.0,
            "dof" => "111",
            "dx" => 0.0,
            "dy" => 0.0,
            "qz" => 0.0
        );
        $this->node[] = array(
            "x" => 120.0,
            "y" => 120.0,
            "dof" => "111",
            "dx" => 0.00,
            "dy" => 0.0,
            "qz" => 0.0
        );
        $this->node[] = array(
            "x" => 120.0,
            "y" => 0.0,
            "dof" => "000",
            "dx" => 0.00,
            "dy" => 0.00,
            "qz" => 0.0
        );
        $this->line[] = array(
            "I" => 0,
            "J" => 1,
            "Mat" => 0,
            "Sec" => 0
        );
        $this->line[] = array(
            "I" => 1,
            "J" => 2,
            "Mat" => 0,
            "Sec" => 1
        );
        $this->line[] = array(
            "I" => 2,
            "J" => 3,
            "Mat" => 0,
            "Sec" => 0
        );
        $this->boundary[] = array(
            "n" => 0,
            "ux" => 0.0,
            "uy" => 0.0,
            "qz" => 0.0
        );
        $this->boundary[] = array(
            "n" => 3,
            "ux" => 0.0,
            "uy" => 0.0,
            "qz" => 0.0
        );
        $this->Mat[] = array(
            "name" => "concrete",
            "E" => 3122.0
        );
        $this->Sec[] = array(
            "name" => "R9x9",
            "A" => 81.0,
            "I" => 546.75
        );
        $this->Sec[] = array(
            "name" => "R9x12",
            "A" => 108.0,
            "I" => 1296.0
        );
        $this->loadcase[0]["nodes"][] = array(
            "n" => 1,
            "Fx" => 1.0,
            "Fy" => 0.0,
            "Mz" => 0.0
            /*
            "F" => array(
                1.0,
                0.0,
                0.0
            )*/
        );
        $this->loadcase[0]["nodes"][] = array(
            "n" => 2,
            "Fx" => 0.0,
            "Fy" => 1.0,
            "Mz" => 0.0
            /*
            "F" => array(
                0.0,
                1.0,
                0.0
            )*/
        );
        $this->loadcase[0]["dists"][] = array(
            "n" => 1,
            "Wx" => 0.0,
            "Wy" => -1.0,
            "Mz" => 0.0
            /*
            "F" => array(
                0.0, -1.0,
                0.0
            )*/
        );
    }

    function Process() // analyse on frame
    
    {
        $ndn = 3;
        $outstr = "<div id='outtext' class='dragobj' style='width:500px;padding:20px'>" . "<h2>Result Table</h2><div style='border:2px solid #333333;padding:10px'>";
        $s = array();
        $f = array();
        $nn = count($this->node);
        $nq = $nn * $ndn;
        $ne = count($this->line);
        $nbw = 6;
        for ($i = 0;$i < $ne;$i++)
        {
            $nb = $ndn * (abs($this->line[$i]["I"] - $this->line[$i]["J"]) + 1.0);
            if ($nbw < $nb) $nbw = $nb;
        }
        //printf($nbw);
        for ($i = 0;$i < $ne;$i++)
        {
            $nmat = $this->line[$i]["Mat"];
            $nsec = $this->line[$i]["Sec"];
            $node1 = $this->node[$this->line[$i]["I"]];
            $node2 = $this->node[$this->line[$i]["J"]];
            $EA = $this->Mat[$nmat]["E"] * $this->Sec[$nsec]["A"];
            $EI = $this->Mat[$nmat]["E"] * $this->Sec[$nsec]["I"];
            $dx = $node2["x"] - $node1["x"];
            $dy = $node2["y"] - $node1["y"];
            //printf(" %3d  %11.4f  %11.4f  %11.4f  %11.4f <br>", $i, $EA, $EI, $dx, $dy);
            //printf(" %3d  %11.4f  %11.4f  %11.4f  %11.4f <br>", $i, $this->Mat[$nmat]["E"], $this->Sec[$nsec]["I"], $dx, $dy);
            $Km = array();
            $this->BmCol2($Km, $EA, $EI, $dx, $dy);
            $node = array(
                $this->line[$i]["I"],
                $this->line[$i]["J"]
            );
            $this->assemble($node, $nbw, $Km, $s);
            //printf("Element %d <br>", $i);
            
        }
        //echo "finished assemble <br/>";
        $cnst = 0;
        for ($i = 0;$i < $nq;$i++)
        {
            if ($cnst < $s[$i * $nbw]) $cnst = $s[$i * $nbw];
        }
        $cnst = 10000 * $cnst;

        $n = count($this->loadcase);
        for ($i = 0;$i < $nq;$i++) $f[$i] = 0;
        for ($i = 0;$i < $n;$i++)
        {
            $this->formf($f, $this->loadcase, $i);
            $this->formelementload($this->node, $this->line, $f, $this->loadcase, $i);
        }
        $this->modifydispbond($this->boundary, $s, $nbw, $cnst, $f);
        $this->bansol($s, $f, $nq, $nbw);
        //echo "solved primary <br/>";
        $outdisp = "<table cellspacing='2px' style='font-size:12px'><tr><th>node#</th><th>x-displ.</th><th>y-displ.</th><th>rotation</th></tr>";
        $outdisp = "<h4>Nodal Displacements</h3>  <div class=\"fixedheadertable\">
        <table border=\"1\" cellpadding=\"5\" cellspacing=\"0\" class=\"whitelinks\">
        <tr><th>node#</th><th>x-displ.</th><th>y-displ.</th><th>Rotation</th></tr>";
        $nn = count($this->node);
        for ($i = 0;$i < $nn;$i++)
        {
            $i1 = $ndn * $i;
            $str = sprintf("<tr><td>%3d</td><td>%01.8f</td><td>%01.8f</td><td>%01.8f</td></tr>", $i, $f[$i1], $f[$i1 + 1], $f[$i1 + 2]);
            $outdisp .= $str;
            $this->node[$i]["dx"] = $f[$i1];
            $this->node[$i]["dy"] = $f[$i1 + 1];
            $this->node[$i]["qz"] = $f[$i1 + 2];
        }
        $outdisp .= "</table></div>";
        //echo "finished Disp <br/>";
        $outstr .= $outdisp;
        $ne = count($this->line);
        
        $outfrc = "<h4>Element Forces</h3>  <div class=\"fixedheadertable\">
        <table border=\"1\" cellpadding=\"5\" cellspacing=\"0\" class=\"whitelinks\">
        <tr><th>member#</th><th>PI</th><th>VI</th><th>RI</th><th>PJ</th><th>VJ</th><th>RJ</th></tr>";
        for ($i = 0;$i < $ne;$i++)
        {
            $nmat = $this->line[$i]["Mat"];
            $nsec = $this->line[$i]["Sec"];
            $node1 = $this->node[$this->line[$i]["I"]];
            $node2 = $this->node[$this->line[$i]["J"]];
            $EA = $this->Mat[$nmat]["E"] * $this->Sec[$nsec]["A"];
            $EI = $this->Mat[$nmat]["E"] * $this->Sec[$nsec]["I"];

            $dx = $node2["x"] - $node1["x"];
            $dy = $node2["y"] - $node1["y"];
            $ELL = sqrt(pow(($dy) , 2) + pow(($dx) , 2));
            //MessageOut(" %3d  %11.4f  %11.4f  %11.4f  %11.4f <br>", $i, $EA, $EI, $dx, $dy);
            //MessageOut(" %3d  %11.4f  %11.4f  %11.4f  %11.4f <br>", $i, $this->Mat[$nmat]["E"], $this->Sec[$nsec]["I"], $dx, $dy);
            $Km = array();
            $this->BmCol1($Km, $EA, $EI, $dx, $dy);
            list($i1, $i2) = array(
                $this->line[$i]["I"] * $ndn,
                $this->line[$i]["J"] * $ndn
            );
            $alambada = $this->lambada($dx, $dy);
            $ed = array();
            for ($j = 0;$j < $ndn;$j++)
            {
                $ed[$j] = $f[$i1 + $j];
                $ed[$j + $ndn] = $f[$i2 + $j];
            }
            for ($j = 0;$j < 6;$j++)
            {
                $edp[$j] = 0;
                for ($k = 0;$k < 6;$k++)
                {
                    $edp[$j] = $edp[$j] + $alambada[$j][$k] * $ed[$k];
                }
            }
            $ed = array();
            for ($j = 0;$j < 6;$j++) $ed[$j] = 0;
            // $ed should member load vector
            $nl = count($this->loadcase);

            for ($j = 0;$j < $nl;$j++)
            {
                $nk = count($this->loadcase[$j]['dists']);
                for ($k = 0;$k < $nk;$k++)
                {
                    $mn = $this->loadcase[$j]['dists'][$k]["n"];
                    if ($i == $mn)
                    {
                        $w = $this->loadcase[$j]['dists'][$k]["Wy"];
                        $ed[0] = 0;
                        $ed[3] = 0;
                        $ed[1] = - $w * $ELL / 2.0;
                        $ed[4] = $ed[1];
                        $ed[2] = - $w * $ELL * $ELL / 12.0;
                        $ed[5] = - $ed[2];
                    }
                }
            }
            for ($j = 0;$j < 6;$j++)
            {
                $ef[$j] = $ed[$j];
                for ($k = 0;$k < 6;$k++)
                {
                    $ef[$j] = $ef[$j] + $Km[$j][$k] * $edp[$k]; // $Km must be local
                    
                }
            }
            $str = sprintf("<tr ><td>%d", $i);
            $outfrc .= $str;
            for ($j = 0;$j < 2;$j++)
            {
                $ii = $ndn * $j;
                $str = sprintf("</td><td>%11.4f</td><td>%11.4f</td><td>%11.4f", $ef[$ii], $ef[$ii + 1], $ef[$ii + 2]);
                $outfrc .= $str;
            }
            $outfrc .= "</td></tr>";
        }
        $outstr .= $outfrc . "</table></div>";
        //echo "finished Element <br/>";
        $this->react = $this->Reaction($this->boundary, $cnst, $f);

        $n = count($this->react);
        
        $outreac = "<h4>Reaction</h3>  <div class=\"fixedheadertable\">
        <table border=\"1\" cellpadding=\"5\" cellspacing=\"0\" class=\"whitelinks\">
        <tr><th>node#</th><th>Rx</th><th>Ry</th><th>Mz</th></tr>";
        for ($i = 0;$i < $n;$i++)
        {
            $rx = $ry = $qz = 0;
            if (isset($this->react[$i]["ux"])) $rx = $this->react[$i]["ux"];
            if (isset($this->react[$i]["uy"])) $ry = $this->react[$i]["uy"];
            if (isset($this->react[$i]["qz"])) $qz = $this->react[$i]["qz"];
            $str = sprintf("<tr><td>%3d</td><td>%01.8f</td><td>%01.8f</td><td>%01.8f</td></tr>", $this->react[$i]["nn"], $rx, $ry, $qz);
            $outreac .= $str;
        };
        $outreac .= "</table></div>";
        //echo "finished React <br/>";
        $outstr .= $outreac . "</div></div>";

        $this->resulthtml = $outstr;
        $this->processfinished = true;
        return $this->resulthtml;
    }
    function Reset()
    {
        $this->node = array();
        $this->line = array();
        $this->boundary = array();
        $this->Mat = array();
        $this->Sec = array();
        $this->nodalmass = array();
        $this->loadcase = array();
    }
    
    function saveCSV($fname){
        $handle = fopen($fname, "w");
        if ($handle)
        {
            $this->saveArr2CSV($handle, $this->node, "node");
            $this->saveArr2CSV($handle, $this->line, "line");
            $this->saveArr2CSV($handle, $this->boundary, "boundary");
            $this->saveArr2CSV($handle, $this->Mat, "Mat");
            $this->saveArr2CSV($handle, $this->Sec, "Sec");
            $this->saveArr2CSV($handle, $this->loadcase[0]["nodes"], "loadcase,0,nodes");
            $this->saveArr2CSV($handle, $this->loadcase[0]["dists"], "loadcase,0,dists");
            fclose($handle);
        }
        else
        {
            // error opening the file.
            
        }
    }
    function readGenCSV($fname){
        $this->Reset();
        Parent::readCSV($fname);
    }
    function readCSVOld($fname)
    {
        $this->Reset();
        $handle = fopen($fname, "r");
        if ($handle)
        {
            while (!feof($handle) && ($line = fgets($handle)) !== false)
            {
                // process the line read.
                if (!empty($line) && $line[0] !== ';')
                {
                    $args = explode(",", $line);
                    switch ($args[0])
                    {
                        case "node":
                            $this->node[] = array(
                                "x" => (float)$args[1],
                                "y" => (float)$args[2],
                                "z" => (float)$args[3],
                                "dof" => $args[4],
                                "dx" => 0.00,
                                "dy" => 0.00,
                                "dz" => 0.0
                            );
                        break;
                        case "line":
                            $this->line[] = array(
                                "I" => (int)$args[1],
                                "J" => (int)$args[2],
                                "Mat" => (int)$args[3],
                                "Sec" => (int)$args[4]
                            );
                        break;
                        case "boundary":
                            $this->boundary[] = array(
                                "n" => (int)$args[1],
                                "ux" => (float)$args[2],
                                "uy" => (float)$args[3],
                                "uz" => (float)$args[4]
                            );
                        break;
                        case "Mat":
                            $this->Mat[] = array(
                                "name" => $args[1],
                                "E" => (float)$args[2],
                                "alpha" => (float)$args[3],
                                "wpv" => (float)$args[4],
                                "mpv" => (float)$args[5]
                            );
                        break;
                        case "Sec":
                            $this->Sec[] = array(
                                "name" => $args[1],
                                "A" => (float)$args[2],
                                "I" => (float)$args[3]
                            );
                        break;
                        case "nodalmass":
                            $this->nodalmass[] = array(
                                "n" => (int)$args[1],
                                "Mx" => (float)$args[2],
                                "My" => (float)$args[3],
                                "Mz" => (float)$args[4]
                            );
                        break;
                        case "loadcase":
                            if ($args[2] === "nodes") $this->loadcase[(int)$args[1]]["nodes"][] = array(
                                "n" => (int)$args[3],
                                "Fx" => (float)$args[4],
                                "Fy" => (float)$args[5],
                                "Fz" => (float)$args[6]
                            );
                            if ($args[2] === "dists") $this->loadcase[(int)$args[1]]["dists"][] = array(
                                "n" => (int)$args[3],
                                "Wx" => (float)$args[4],
                                "Wy" => (float)$args[5],
                                "Wz" => (float)$args[6]
                            );
                            break;
                        }
                    }

            }
            fclose($handle);
        }
        else
        {
            // error opening the file.
            
        }
    }
    function ncolbeam()
    {
        $ne = count($this->line);
        $nbeam = 0;
        $ncol = 0;
        for ($i = 0;$i < $ne;$i++)
        {
            $node1 = $this->node[$this->line[$i]["I"]];
            $node2 = $this->node[$this->line[$i]["J"]];
            $dx = $node2["x"] - $node1["x"];
            $dy = $node2["y"] - $node1["y"];
            if ($dx == 0) $ncol++;
            else if ($dy == 0) $nbeam++;
        }
        return array(
            $ncol,
            $nbeam
        );
    }
    function maxreaction()
    {
        $maxr = 0;
        $nn = 0;
        $n = count($this->boundary);
        for ($i = 0;$i < $n;$i++)
        {
            $r = $this->react[$i]["uy"];
            if (abs($maxr) < abs($r))
            {
                $maxr = $r;
                $nn = $i;
            }
        }
        return array(
            $nn,
            $maxr
        );
    }
    function maxdisp()
    {
        $maxr = 0;
        $nn = 0;
        $n = count($this->node);
        for ($i = 0;$i < $n;$i++)
        {
            $dx = $this->node[$i]["dx"];
            $dy = $this->node[$i]["dy"];
            $dl = sqrt($dx * $dx + $dy * $dy);
            if ($maxr < $dl)
            {
                $maxr = $dl;
                $nn = $i;
            }
        }
        return array(
            $nn,
            $maxr
        );
    }
}
//echo " pframe";

?>

<?php
/*
$pframe = new Frame2D();
$pframe->Init();
//echo " init";
echo $pframe->Process();
$pframe->saveCSV("frame2d.csv");
$pframe->readGenCSV("frame2d.csv");
echo $pframe->Process();
*/
?>

                                                