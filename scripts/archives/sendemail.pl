#!/usr/bin/perl -w

###########################################################################
##
## send_mail.pl Copyright (c) 1998-2004, Peter M. Kekenes Huskey, 
## All Rights Reserved.
##
## This script may be modified without the express permission of the    
## author.  However this notice must remain intact.                     
##
###########################################################################	

#
# peter m kekenes-huskey
#
# 2003
#

use strict;
use Net::SMTP;
use CGI; 
print header();

print "THIS IS A TEST\n";
die "\n";

#
#  Main
#
{
  # Globals
  my $user = "tod";
  my $email_address = "${user}\@its.caltech.edu";
  my $SMTP_ServerName = "mail-fwd.sbc-webhosting.com";
	
  SendEmail($email_address, $SMTP_ServerName);

}


########################################################
# Establish connection with SMTP server and send email #
########################################################
sub SendEmail {
  my ($email_address, $SMTP_ServerName) = @_;

  print "Sending email to $email_address\n";
  
  # Connect to the server
  my $smtp = Net::SMTP->new($SMTP_ServerName);
  die "Couldn't connect to server" unless $smtp;
  
  #define recipients
  #my $sender_address = $email_address;
  my $sender_address = "csupport\@landtemp.com";
  my $recipient_address = $email_address;
  $smtp->mail($sender_address);
  $smtp->to($recipient_address);
  
  
  # Start the mail
  $smtp->data();
  
  # Send the header.
  $smtp->datasend("To: $recipient_address\n");
  $smtp->datasend("From: $sender_address\n");
  $smtp->datasend("Subject: test\n");
  $smtp->datasend("\n");
  
  # Send the message
  $smtp->datasend("Can you do stat mech on tomorrow or friday?\n\n");
  
  # Send the termination string
  $smtp->dataend();
  
  $smtp->quit();

}

