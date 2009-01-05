#!/usr/bin/env ruby

$: << File.join(File.expand_path(File.dirname(__FILE__)), %w(.. lib))

require 'admix'
require 'admixweb/prettify'

require 'rubygems'
require 'sinatra'

module Admix

  get '/' do
    haml :index
  end

  post '/admix' do |*args|
    begin
      loc = params[:loc][:tempfile].read
      ped = params[:ped][:tempfile].read
      result = Wrapper.call :loc => loc, :ped => ped
      case params[:type]
      when "html"
        content_type 'text/html'
        pretty = AdmixWeb::prettify result
        haml :results, :locals => { :pretty => pretty }
      when "text"
        content_type 'text/plain'
        result
      else
        throw :halt, [500, "Invalid form"]
      end
    rescue AdmixError => e
      throw :halt, [500, "Error: #{e.message}"]
    end
  end

  get "/style.css" do
    content_type "text/css"

    <<EOCSS
form dd { margin: 1ex 0 1ex 15ex }
form dt { float: left;
          text-align: right;
          width: 14ex;
          margin: 0; }
EOCSS
  end

  template :index do <<EOHAML
!!!
%title Admix
%link{ :rel => "stylesheet", :type => "text/css", :href => "style.css" }
%form{ :action => "/admix", :method => "post", :enctype => "multipart/form-data"}
  %dl
    %dt
      %label{ :for => :loc } Locus file
    %dd
      %input{ :name => :loc, :id => :loc, :type => :file }
    %dt
      %label{ :for => :ped } Pedigree file
    %dd
      %input{ :name => :ped, :id => :ped, :type => :file }
    %dd
      %label
        %input{ :type => :radio, :name => :type, :value => :html }
        HTML
      %label
        %input{ :type => :radio, :name => :type, :value => :text, :checked => true }
        Plain text
    %dd
      %button Submit
EOHAML
  end

  template :results do <<EOHAML
!!!
%title Admix Results
%table= pretty
EOHAML
  end

  configure do
    set_option :haml, :format => :html5
  end

end
