function header(){
  new Ajax.Updater("header","/_header",{ method:"GET" });
}

function footer(){
  new Ajax.Updater("footer","/_footer",{ method:"GET" });
}
