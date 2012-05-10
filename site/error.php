<?php
	$title="Error";
	include('includes/page-top.php');
?>

<div id="main">
<h1>Error</h1>
<?php
switch(getenv("REDIRECT_STATUS"))
{
	# "400 Bad Request"
	case 400:
		$error = "Bad Server Request (Error 400)";
		break;

	# "401 Authorization Required"
	case 401:
		$error = "Authorization Required (Error 401)";
		$explain = "This section requires a password or is otherwise protected. If you feel you have received this message in error, please contact your friendly neighbourhood admin for access.";
		break;

# "403 Forbidden"
	case 403:
		$error = "Forbidden Access (Error 403)" ;
		$explain = "This section requires a password or is otherwise protected. You may not pass.";
		break;

	# "404 Not Found"
	case 404:
		$error = "Not found (Error 404)";
		$explain = "We've recently switched to a new site and some things may be missing or broken. Please let us know what you were trying to access, and we'll find it for you! Apologies for the inconvenience.";
		break;

	# "500 Internal Server Error"
	case 500:
		$error = "Internal Server Error (Error 500)";
		$explain = "Please check the address and try again.";
		break;
}
print "<h1>$error</h1>\n";
print "<p>$explain</p>\n";
?>

</div>

<?php include('includes/page-bottom.php'); ?>
