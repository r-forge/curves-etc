
<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->

<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='https://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"https://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="https://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">

  <head>
	<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
	<title><?php echo $group_name; ?></title>
	<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
  </head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="/"><img src="<?php echo $themeroot; ?>/images/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
</table>


<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->

<?php if ($handle=fopen('https://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>

<!-- end of project description -->
<p>
This is the Subversion repository of <strong>several</strong> <a href="https://www.r-project.org">R</a> packages, all related to <em>curve estimation</em> aka curve <em>"fitting"</em>.
</p>

<ul>
<li> The individual packages' CRAN pages are
<p>
<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">
<thead>
<tr><th scope="col" class="left">Name</th><th scope="col" class="left">CRAN page</th><th scope="col" class="left">Short Description</th></tr>
</thead>
<tbody>
<tr><td>cobs</td><td><a href="https://cran.r-project.org/package=cobs">https://cran.r-project.org/package=cobs</a></td><td>Constrained B-Splines</td></tr>
<tr><td>cobs99</td><td><a href="https://cran.r-project.org/package=cobs99">https://cran.r-project.org/package=cobs99</a></td><td>COBS - old 1999 version</td></tr>
<tr><td>lokern</td><td><a href="https://cran.r-project.org/package=lokern">https://cran.r-project.org/package=lokern</a></td><td>Kernel Regression Smoothing, Local/Global Bandwidth</td></tr>
<tr><td>lpridge</td><td><a href="https://cran.r-project.org/package=lpridge">https://cran.r-project.org/package=lpridge</a></td><td>Local Polynomial (Ridge) Regression</td></tr>
<tr><td>nor1mix</td><td><a href="https://cran.r-project.org/package=nor1mix">https://cran.r-project.org/package=nor1mix</a></td><td>Normal (1-d) Mixture Models</td></tr>
<tr><td>plugdensity</td><td><a href="https://cran.r-project.org/package=plugdensity">https://cran.r-project.org/package=plugdensity</a></td><td>Hermann`s Plugin Density Estimate</td></tr>
</tbody>
</table></p>
</li>

<li> The <strong>R-forge project summary page</strong> you can
  find <a href="https://<?php echo $domain; ?>/projects/<?php echo
		$group_name; ?>/"><strong>here</strong></a>.
</li>
</ul>

</body>
</html>
