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
      result = Wrapper.call :loc => params[:loc], :ped => params[:ped]
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

  template :index do <<EOHAML
!!!
%title Admix
%form{ :action => "/admix", :method => "post", :enctype => "multipart/form-data"}
  %label
    Locus file
    %input{ :name => :loc, :type => "file" }
  %label
    Pedigree file
    %input{ :name => :ped, :type => "file" }
  %label
    %input{ :type => :radio, :name => :type, :value => :html }
    HTML
  %label
    %input{ :type => :radio, :name => :type, :value => :text, :checked => true }
    Plain text
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
