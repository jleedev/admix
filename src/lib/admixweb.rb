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
        AdmixWeb::prettify result
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
%form{ :action => "/admix", :method => "post"}
  %label
    Locus file
    %textarea{ :name => :loc }
  %label
    Pedigree file
    %textarea{ :name => :ped }
  %label
    %input{ :type => :radio, :name => :type, :value => :html }
    HTML
  %label
    %input{ :type => :radio, :name => :type, :value => :text, :checked => true }
    Plain text
  %button Submit
EOHAML
  end

  configure do
    set_option :haml, :format => :html5
  end

end
